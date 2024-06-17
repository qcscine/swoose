/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include <Core/Log.h>
#include <Molassembler/Graph.h>
#include <Molassembler/Interpret.h>
#include <Molassembler/Molecule.h>
#include <Swoose/MMParametrization/MMParametrizationSettings.h>
#include <Swoose/MMParametrization/MolecularSystemPartitioner.h>
#include <Swoose/MMParametrization/ParametrizationData.h>
#include <Swoose/MMParametrization/ParametrizationUtils/AtomicChargesAssembler.h>
#include <Swoose/MMParametrization/ParametrizationUtils/ConnectivityGenerator.h>
#include <Swoose/MMParametrization/ParametrizationUtils/FullHessianAssembler.h>
#include <Swoose/MMParametrization/ParametrizationUtils/SuperfluousFragmentIdentifier.h>
#include <Swoose/MolecularMechanics/SFAM/SfamAtomTypeIdentifier.h>
#include <Swoose/MolecularMechanics/Topology/IndexedStructuralTopologyCreator.h>
#include <Utils/Bonds/BondDetector.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>
#include <numeric>

using namespace testing;
namespace Scine {
using namespace MMParametrization;
namespace Tests {

/**
 * @class AParametrizationOfLargeSystemsTest ParametrizationOfLargeSystemsTest.cpp
 * @brief Tests the automated parametrization for large molecular systems.
 * @test
 */
class AParametrizationOfLargeSystemsTest : public Test {
 public:
  void SetUp() override {
    boost::filesystem::remove_all(connFile);
  }
  std::string connFile = "Connectivity.dat";
  ParametrizationData data;
  Core::Log silentLogger = Core::Log::silent();
};

TEST_F(AParametrizationOfLargeSystemsTest, FragmentOfDermcidinIsFragmentedCorrectly) {
  // This test does not act on the full dermcidin system, just on a previously generated fragment of it
  // due to the fact that in debug mode this test would take too much time.
  Utils::AtomCollection largeStructure = Utils::ChemicalFileHandler::read(dermcidin_fragment_xyz_file).first;
  int numberOfAtoms = largeStructure.size();

  // Initialize and modify settings
  auto settings = std::make_shared<MMParametrizationSettings>();
  settings->modifyInt(SwooseUtilities::SettingsNames::numberAtomsThreshold, 120);
  settings->modifyDouble(SwooseUtilities::SettingsNames::subsystemRadius, 5.0); // Rather small radius for this test

  // Move structure to data object
  data.fullStructure = std::move(largeStructure);
  // Set number of atoms in data object
  data.numberOfAtoms = data.fullStructure.size();

  // Generate connectivity and update the listsOfNeighbors member of data_
  ConnectivityGenerator connectivityGenerator(data, settings, silentLogger);
  connectivityGenerator.generateInitialListsOfNeighbors();

  data.unpairedElectrons[88] = 1; // Build in intentional mistake
  // Partition system into subsystems
  MolecularSystemPartitioner partitioner(data, settings, silentLogger);
  std::string exceptionString = "";
  try {
    partitioner.divideSystem();
  }
  catch (const std::runtime_error& e) {
    exceptionString = e.what();
  }

  ASSERT_STREQ(exceptionString.c_str(), "Error while trying to generate fragment with index 0.");

  data.unpairedElectrons[88] = 0; // Correct intentional mistake
  partitioner.divideSystem();     // Now start the correct fragmentation run

  // Number of fragments should be the same as the number of atoms in the molecular system
  ASSERT_THAT(data.vectorOfStructures.size(), Eq(numberOfAtoms));
  for (const auto& s : data.vectorOfStructures) {
    ASSERT_THAT(s->size(), Gt(19));  // All fragments have at least 20 atoms
    ASSERT_THAT(s->size(), Lt(100)); // All fragments have less than 100 atoms

    int totalNumberOfElectrons = 0;
    for (int i = 0; i < s->size(); ++i) {
      totalNumberOfElectrons += Utils::ElementInfo::Z(s->getElement(i));
    }
    // All fragments should have an even number of electrons for this molecular system
    ASSERT_TRUE(totalNumberOfElectrons % 2 == 0);

    // Get information about the individual molecules in a fragment with molassembler
    auto bondOrders = Utils::BondDetector::detectBonds(*s);
    auto result = Molassembler::Interpret::graphs(*s, bondOrders);

    // There should not be a molecule in a fragment that has 4 atoms or fewer.
    for (const auto& g : result.graphs) {
      auto nAtoms = g.V();
      ASSERT_TRUE(nAtoms > 4);
    }
  }

  // Assert that the constrained atoms were correctly set
  for (int k = 0; k < int(data.constrainedAtoms.size()); ++k) {
    auto constrainedAtoms = data.constrainedAtoms[k];
    ASSERT_TRUE(constrainedAtoms.empty() || constrainedAtoms.size() > 2);

    // Check some specific cases
    if (k == 0) {
      ASSERT_THAT(constrainedAtoms.size(), Eq(10));
      ASSERT_THAT(constrainedAtoms[0], Eq(1));
      ASSERT_THAT(constrainedAtoms[1], Eq(54));
      ASSERT_THAT(constrainedAtoms[2], Eq(5));
      ASSERT_THAT(constrainedAtoms[3], Eq(56));
    }
    else if (k == 58 || k == 104) {
      ASSERT_TRUE(constrainedAtoms.empty());
    }
    else if (k == 115) {
      ASSERT_THAT(constrainedAtoms.size(), Eq(4));
      ASSERT_THAT(constrainedAtoms[0], Eq(7));
      ASSERT_THAT(constrainedAtoms[1], Eq(26));
      ASSERT_THAT(constrainedAtoms[2], Eq(24));
      ASSERT_THAT(constrainedAtoms[3], Eq(27));
    }
  }

  // Test for one specific system that the fragment is correct
  Utils::AtomCollection testFragment = *data.vectorOfStructures[46];
  ASSERT_THAT(testFragment.size(), Eq(37));
  Eigen::RowVector3d positionOfOneSaturatingHydrogen = testFragment.getPosition(32);
  ASSERT_THAT(positionOfOneSaturatingHydrogen(0), DoubleNear(5.88945, 1e-3));
  ASSERT_THAT(positionOfOneSaturatingHydrogen(1), DoubleNear(14.8917, 1e-3));
  ASSERT_THAT(positionOfOneSaturatingHydrogen(2), DoubleNear(20.6470, 1e-3));

  // Sanity check of that test fragment
  for (int k = 0; k < testFragment.size(); ++k) {
    for (int l = 0; l < k; ++l) {
      auto distance = (testFragment.getPosition(k) - testFragment.getPosition(l)).norm();
      ASSERT_TRUE(distance > 1.8); // No atoms closer than 1.8 bohr
    }
  }

  // Testing the superfluous fragment identifier
  SuperfluousFragmentIdentifier::identifySuperfluousFragments(data, silentLogger, nullptr);
  ASSERT_THAT(data.superfluousFragments.size(), Eq(17));
}

TEST_F(AParametrizationOfLargeSystemsTest, TopologyIsGeneratedCorrectlyForDermcidin) {
  Utils::AtomCollection largeStructure = Utils::ChemicalFileHandler::read(dermcidin_xyz_file).first;
  auto settings = std::make_shared<MMParametrizationSettings>();

  // Move structure to data object
  data.fullStructure = std::move(largeStructure);
  // Set number of atoms in data object
  data.numberOfAtoms = data.fullStructure.size();

  // Generate connectivity and update the listsOfNeighbors member of data_
  ConnectivityGenerator connectivityGenerator(data, settings, silentLogger);
  connectivityGenerator.generateInitialListsOfNeighbors();

  // Generate topology
  MolecularMechanics::IndexedStructuralTopologyCreator topologyCreator(data.listsOfNeighbors);
  data.topology = topologyCreator.calculateIndexedStructuralTopology();
  topologyCreator.addHydrogenBondsToIndexedStructuralTopology(data.topology, data.fullStructure);

  ASSERT_THAT(data.topology.getBondContainer().size(), Eq(697));
  ASSERT_THAT(data.topology.getAngleContainer().size(), Eq(1263));
  ASSERT_THAT(data.topology.getDihedralContainer().size(), Eq(1832));
  ASSERT_THAT(data.topology.getImproperDihedralContainer().size(), Eq(116));
  ASSERT_THAT(data.topology.getExcludedNonBondedContainer().size(), Eq(1960));
  ASSERT_THAT(data.topology.getScaledNonBondedContainer().size(), Eq(1827));
  ASSERT_THAT(data.topology.getHydrogenBondContainer().size(), Eq(745));
}

TEST_F(AParametrizationOfLargeSystemsTest, FullHessianAssemblerWorksCorrectly) {
  // Try this just on a small system (e.g., alanine)
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;
  auto settings = std::make_shared<MMParametrizationSettings>();

  // Move structure to data object
  data.fullStructure = std::move(structure);
  // Set number of atoms in data object
  data.numberOfAtoms = data.fullStructure.size();

  // Generate connectivity and update the listsOfNeighbors member of data
  ConnectivityGenerator connectivityGenerator(data, settings, silentLogger);
  connectivityGenerator.generateInitialListsOfNeighbors();

  // Generate topology
  MolecularMechanics::IndexedStructuralTopologyCreator topologyCreator(data.listsOfNeighbors);
  data.topology = topologyCreator.calculateIndexedStructuralTopology();

  // Next several lines are used later for testing two specific subblocks for correctness
  Eigen::MatrixXd subblockOriginalOne(3, 3);
  Eigen::MatrixXd subblockOriginalTwo(3, 3);
  Eigen::MatrixXd subblockOriginalThree(3, 3);
  Eigen::MatrixXd subblockOriginalFour(3, 3);
  const int index5 = 5;
  const int index7 = 7;
  const int index8 = 8;
  const int index9 = 9;
  const int index10 = 10;
  const int index11 = 11;
  const int index12 = 12;

  // Random Hessians for all subsystems
  data.vectorOfStructures.reserve(data.numberOfAtoms);
  data.vectorOfHessians.reserve(data.numberOfAtoms);
  for (int i = 0; i < data.numberOfAtoms; ++i) {
    // Each subsystem is going to be the whole alanine molecule
    data.vectorOfStructures.push_back(std::make_unique<Utils::AtomCollection>(data.fullStructure));
    if (i > 1) {
      Utils::HessianMatrix hessian = Eigen::MatrixXd::Random(3 * data.numberOfAtoms, 3 * data.numberOfAtoms);

      if (i == 5) {
        subblockOriginalOne = hessian.block<3, 3>(3 * index5, 3 * index7);
        subblockOriginalTwo = hessian.block<3, 3>(3 * index7, 3 * index8);
      }
      else if (i == 9) {
        subblockOriginalThree = hessian.block<3, 3>(3 * index9, 3 * index11);
        subblockOriginalFour = hessian.block<3, 3>(3 * index10, 3 * index12);
      }

      data.vectorOfHessians.push_back(std::make_unique<Utils::HessianMatrix>(hessian));
    }
    else {
      // Test whether the FullHessianAssembler can handle that the Hessian for atoms 0 and 1 is not available
      data.vectorOfHessians.push_back(nullptr);
    }
    // Populate the atom index mapping vector
    std::vector<int> indicesOfAtomsInFragment(data.numberOfAtoms);
    std::iota(indicesOfAtomsInFragment.begin(), indicesOfAtomsInFragment.end(), 0);
    data.atomIndexMapping.push_back(indicesOfAtomsInFragment);
  }

  FullHessianAssembler hessianAssembler(data, silentLogger);
  hessianAssembler.assembleFullHessian();

  double densityInPercent = 100.0 * data.fullHessian.nonZeros() / data.fullHessian.size();
  ASSERT_TRUE(densityInPercent < 40);

  // Block 5 - 7 should come from fragment 5 (due to bond: 5 - 7)
  Eigen::MatrixXd subblockFromSparseHessianOne = data.fullHessian.block(3 * index5, 3 * index7, 3, 3);

  // Block 7 - 8 should come from fragment 5 (due to angle: 7 - 5 - 8)
  Eigen::MatrixXd subblockFromSparseHessianTwo = data.fullHessian.block(3 * index7, 3 * index8, 3, 3);

  // Block 9 - 11 should come from fragment 9 (due to bond: 9 - 11)
  Eigen::MatrixXd subblockFromSparseHessianThree = data.fullHessian.block(3 * index9, 3 * index11, 3, 3);

  // Block 10 - 12 should come from fragment 9 (due to dihedral: 10 - 9 - 11 - 12)
  Eigen::MatrixXd subblockFromSparseHessianFour = data.fullHessian.block(3 * index10, 3 * index12, 3, 3);

  for (int k = 0; k < 3; ++k) {
    for (int l = 0; l < 3; ++l) {
      ASSERT_THAT(subblockOriginalOne(k, l), DoubleNear(subblockFromSparseHessianOne(k, l), 1e-6));
      ASSERT_THAT(subblockOriginalTwo(k, l), DoubleNear(subblockFromSparseHessianTwo(k, l), 1e-6));
      ASSERT_THAT(subblockOriginalThree(k, l), DoubleNear(subblockFromSparseHessianThree(k, l), 1e-6));
      ASSERT_THAT(subblockOriginalFour(k, l), DoubleNear(subblockFromSparseHessianFour(k, l), 1e-6));
    }
  }
}

TEST_F(AParametrizationOfLargeSystemsTest, AtomicChargesAssemblerWorksCorrectly) {
  Utils::AtomCollection largeStructure = Utils::ChemicalFileHandler::read(dermcidin_xyz_file).first;
  auto settings = std::make_shared<MMParametrizationSettings>();
  data.fullStructure = largeStructure;
  data.numberOfAtoms = largeStructure.size();
  data.vectorOfStructures.resize(data.numberOfAtoms);

  // Generate connectivity and update the listsOfNeighbors member of data
  ConnectivityGenerator connectivityGenerator(data, settings, silentLogger);
  connectivityGenerator.generateInitialListsOfNeighbors();

  // Generate topology and atom types
  MolecularMechanics::IndexedStructuralTopologyCreator topologyCreator(data.listsOfNeighbors);
  data.topology = topologyCreator.calculateIndexedStructuralTopology();
  topologyCreator.addHydrogenBondsToIndexedStructuralTopology(data.topology, data.fullStructure);
  MolecularMechanics::SfamAtomTypeIdentifier atomTypeIdentifier(data.numberOfAtoms, data.fullStructure.getElements(),
                                                                data.listsOfNeighbors);
  MolecularMechanics::SfamAtomTypeLevel atl =
      MolecularMechanics::SfamAtomTypeIdentifier::generateSfamAtomTypeLevelFromString("high");
  data.atomTypes = atomTypeIdentifier.getAtomTypes(atl);

  // Generate a mock vector of atomic charges for each fragment
  data.atomicChargesForEachFragment.reserve(largeStructure.size());
  for (int i = 0; i < largeStructure.size(); ++i) {
    std::vector<double> atomicCharges;
    for (int k = 0; k < 50; ++k) {
      if (k % 2 == 0)
        atomicCharges.push_back(0.9);
      else
        atomicCharges.push_back(-0.1);
    }
    data.atomicChargesForEachFragment.push_back(atomicCharges);

    // Generate a mock index map
    // The following index map will result in all atoms with index larger 49 will have a charge of -0.1 and
    // the others will have alternating charges of 0.9 and -0.1
    std::vector<int> mockIndexMap(49);
    std::iota(mockIndexMap.begin(), mockIndexMap.end(), 0);
    if (std::find(mockIndexMap.begin(), mockIndexMap.end(), i) == mockIndexMap.end())
      mockIndexMap.push_back(i);
    else
      mockIndexMap.push_back(50);
    data.atomIndexMapping.push_back(mockIndexMap);
  }

  AtomicChargesAssembler::assembleAtomicCharges(data, silentLogger);
  double sumOfAtomicCharges = std::accumulate(data.atomicCharges.begin(), data.atomicCharges.end(), 0.0);
  ASSERT_THAT(sumOfAtomicCharges, DoubleNear(-44.7, 1e-6));
  ASSERT_THAT(data.atomicCharges[2], DoubleNear(0.9, 1e-6));
  ASSERT_THAT(data.atomicCharges[11], DoubleNear(-0.1, 1e-6));
  ASSERT_THAT(data.atomicCharges[117], DoubleNear(-0.1, 1e-6));
  ASSERT_THAT(data.atomicCharges[118], DoubleNear(-0.1, 1e-6));

  // Renormalize the charges
  AtomicChargesAssembler::renormalizeAtomicCharges(data);
  sumOfAtomicCharges = std::accumulate(data.atomicCharges.begin(), data.atomicCharges.end(), 0.0);
  ASSERT_THAT(sumOfAtomicCharges, DoubleNear(0.0, 1e-6));
  ASSERT_THAT(data.atomicCharges[2], DoubleNear(0.9641, 1e-4));
  ASSERT_THAT(data.atomicCharges[11], DoubleNear(-0.0358, 1e-4));
  ASSERT_THAT(data.atomicCharges[117], DoubleNear(-0.0358, 1e-4));
  ASSERT_THAT(data.atomicCharges[118], DoubleNear(-0.0358, 1e-4));

  auto copyOfAtomicChargesVector = data.atomicChargesForEachFragment;

  data.atomicChargesForEachFragment[500].clear();
  // This should throw since atom 500 will not be in the fragment corresponding to one of its covalent neighbors
  EXPECT_THROW(AtomicChargesAssembler::assembleAtomicCharges(data, silentLogger), std::runtime_error);

  data.atomicChargesForEachFragment = copyOfAtomicChargesVector;

  data.atomicChargesForEachFragment[10].clear();
  // Should be compensated and no throw is expected
  AtomicChargesAssembler::assembleAtomicCharges(data, silentLogger);

  for (const auto& neighbor : data.listsOfNeighbors[10]) {
    data.atomicChargesForEachFragment[neighbor].clear();
    for (const auto& neighborOfNeighbor : data.listsOfNeighbors[neighbor])
      data.atomicChargesForEachFragment[neighborOfNeighbor].clear();
  }

  // This should throw again, since the absence of all those charge vectors cannot be compensated anymore
  EXPECT_THROW(AtomicChargesAssembler::assembleAtomicCharges(data, silentLogger), std::runtime_error);
}

} // namespace Tests
} // namespace Scine
