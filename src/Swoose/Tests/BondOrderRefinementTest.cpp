/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include <Core/Log.h>
#include <Swoose/MMParametrization/MMParametrizationSettings.h>
#include <Swoose/MMParametrization/ParametrizationData.h>
#include <Swoose/MMParametrization/ParametrizationUtils/ConnectivityGenerator.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>
#include <numeric>

using namespace testing;
namespace Scine {
using namespace MMParametrization;
namespace Tests {

/**
 * @class BondOrderRefinementTest BondOrderRefinementTest.cpp
 * @brief Tests the refinement of bond orders by quantum-chemical data
 * @test
 */
class BondOrderRefinementTest : public Test {
 public:
  void SetUp() override {
    if (boost::filesystem::exists(connFile))
      boost::filesystem::remove_all(connFile);
  }
  std::string connFile = "Connectivity.dat";
};

TEST_F(BondOrderRefinementTest, ListsOfNeighborsAreRefinedCorrectly) {
  // Create silent logger
  auto silentLogger = Core::Log::silent();
  // Initiate parametrization data and settings objects
  ParametrizationData data;
  auto settings = std::make_shared<MMParametrizationSettings>();

  // Read in full structure
  data.fullStructure = Utils::ChemicalFileHandler::read(dermcidin_fragment_xyz_file).first;
  data.numberOfAtoms = data.fullStructure.size();

  // Fill vector of structures
  data.vectorOfStructures.resize(data.numberOfAtoms);
  for (int i = 0; i < data.numberOfAtoms; ++i) {
    data.vectorOfStructures[i] = std::make_unique<Utils::AtomCollection>(data.fullStructure);
    // Populate the atom index mapping vector
    std::vector<int> indicesOfAtomsInFragment(data.numberOfAtoms);
    std::iota(indicesOfAtomsInFragment.begin(), indicesOfAtomsInFragment.end(), 0);
    data.atomIndexMapping.push_back(indicesOfAtomsInFragment);
  }

  // Generate initial lists of neighbors
  ConnectivityGenerator connectivityGeneratorNormalThreshold(data, settings, silentLogger);
  connectivityGeneratorNormalThreshold.generateInitialListsOfNeighbors();
  bool foundBondBetweenSixAndSeven =
      (std::find(data.listsOfNeighbors[7].begin(), data.listsOfNeighbors[7].end(), 6) != data.listsOfNeighbors[7].end());
  bool foundBondBetweenFourtyAndNinety =
      (std::find(data.listsOfNeighbors[40].begin(), data.listsOfNeighbors[40].end(), 90) != data.listsOfNeighbors[40].end());
  ASSERT_TRUE(foundBondBetweenSixAndSeven);
  ASSERT_FALSE(foundBondBetweenFourtyAndNinety);
  // Make a copy of the initial lists of neighbors for later
  std::vector<std::list<int>> copyOfInitialListsOfNeighbors = data.listsOfNeighbors;

  // Fill vector of bond order collections with only null pointers
  data.vectorOfBondOrderCollections.resize(data.numberOfAtoms);
  for (int i = 0; i < data.numberOfAtoms; ++i) {
    data.vectorOfBondOrderCollections[i] = nullptr;
  }

  // Test out new connectivity generator instance with modified settings
  settings->modifyDouble(SwooseUtilities::SettingsNames::bondOrderThreshold, 0.8); // unusually large threshold
  ConnectivityGenerator connectivityGeneratorLargeThreshold(data, settings, silentLogger);

  // Generate new vector of bond order collections
  for (int i = 0; i < data.numberOfAtoms; ++i) {
    Eigen::MatrixXd bondOrderMatrix = Eigen::MatrixXd::Random(data.numberOfAtoms, data.numberOfAtoms);
    if (i == 6) {
      bondOrderMatrix(6, 7) = 0.7;
      bondOrderMatrix(7, 6) = 0.7;
    }
    if ((i == 40) || (i == 90)) {
      bondOrderMatrix(40, 90) = 0.9;
      bondOrderMatrix(90, 40) = 0.9;
    }
    Utils::BondOrderCollection bondOrders;
    bondOrders.setMatrix(bondOrderMatrix.sparseView());
    data.vectorOfBondOrderCollections[i] = std::make_unique<Utils::BondOrderCollection>(std::move(bondOrders));
    if ((i == 7) || (i == 32) || (i == 70))
      data.vectorOfBondOrderCollections[i] = nullptr;
  }

  data.listsOfNeighbors = copyOfInitialListsOfNeighbors;
  connectivityGeneratorLargeThreshold.refineListsOfNeighbors();
  foundBondBetweenSixAndSeven =
      (std::find(data.listsOfNeighbors[7].begin(), data.listsOfNeighbors[7].end(), 6) != data.listsOfNeighbors[7].end());
  foundBondBetweenFourtyAndNinety = (std::find(data.listsOfNeighbors[40].begin(), data.listsOfNeighbors[40].end(), 90) !=
                                     data.listsOfNeighbors[40].end());
  ASSERT_FALSE(foundBondBetweenSixAndSeven);
  ASSERT_TRUE(foundBondBetweenFourtyAndNinety);
}

} // namespace Tests
} // namespace Scine
