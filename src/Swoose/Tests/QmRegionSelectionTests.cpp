/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Swoose/MolecularMechanics/SFAM/SfamMolecularMechanicsCalculator.h>
#include <Swoose/QMMM/QmRegionSelection/QmRegionCandidateGenerator.h>
#include <Swoose/QMMM/QmRegionSelection/QmRegionSelector.h>
#include <Swoose/QMMM/QmRegionSelection/QmRegionSelectorSettings.h>
#include <Swoose/QMMM/QmRegionSelection/QmmmModelAnalyzer.h>
#include <Swoose/QMMM/QmRegionSelection/QmmmReferenceDataManager.h>
#include <Swoose/Utilities/CalculatorOptions.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/FragmentationHelper.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/Yaml.h>
#include <gmock/gmock.h>
#include <yaml-cpp/yaml.h>
#include <numeric>

using namespace testing;
namespace Scine {
using namespace Qmmm;
namespace Tests {

/**
 * @class QmRegionSelectionTests QmRegionSelectionTests.cpp
 * @brief Tests the automated QM region selection.
 * @test
 */
class QmRegionSelectionTests : public Test {
 public:
  Core::Log silentLogger = Core::Log::silent();

  void SetUp() override {
    srand(42);
  }
};

TEST_F(QmRegionSelectionTests, SingleQmRegionIsCorrectlyConstructed) {
  QmRegionSelector qmRegionSelector;
  qmRegionSelector.setLog(silentLogger);
  qmRegionSelector.settings().modifyIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms, {302});
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::initialRadiusForQmRegionSelection, 5.0);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMinSize, 60);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMaxSize, 70);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionRefMaxSize, 120);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::numAttemptsPerRadius, 3);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::maxNumRefModels, 3);
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::cuttingProbability, 0.9);
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, insulin_connectivity_file);
  /*
   * The following is to induce the correct exception below.
   * However, this might soon change as the other reference
   * data modes become available in the implementation.
   */
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode,
                                           SwooseUtilities::OptionNames::readFromFilesMode);

  EXPECT_THROW(qmRegionSelector.getQmRegionStructure(), QmRegionHasNotBeenSelectedException);
  EXPECT_THROW(qmRegionSelector.getQmRegionIndices(), QmRegionHasNotBeenSelectedException);
  EXPECT_THROW(qmRegionSelector.getQmRegionChargeAndMultiplicity(), QmRegionHasNotBeenSelectedException);

  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(insulin_xyz_file).first;
  std::string exceptionString = "";
  try {
    qmRegionSelector.generateQmRegion(structure);
  }
  catch (const std::runtime_error& e) {
    exceptionString = e.what();
  }
  ASSERT_STREQ(exceptionString.c_str(),
               "Currently, only database mode and direct mode are implemented for calculating QM/MM model forces.");

  // Again with correct cutting probability
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::cuttingProbability, 1.0);
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::initialRadiusForQmRegionSelection, 7.5);
  qmRegionSelector.generateQmRegion(structure);

  auto qmRegion = qmRegionSelector.getQmRegionStructure();
  auto qmRegionIndices = qmRegionSelector.getQmRegionIndices();
  auto qmRegionInfo = qmRegionSelector.getQmRegionChargeAndMultiplicity();

  ASSERT_THAT(qmRegion.size(), Eq(59));
  ASSERT_THAT(qmRegionIndices.size(), Eq(59));
  ASSERT_THAT(qmRegionInfo.first, Eq(0));
  ASSERT_THAT(qmRegionInfo.second, Eq(1));

  std::array<int, 3> atomsInside = {300, 302, 322};
  std::array<int, 5> atomsOutside = {50, 81, 92, 102, 186};
  for (int idx : atomsInside) {
    bool isInside = std::find(qmRegionIndices.begin(), qmRegionIndices.end(), idx) != qmRegionIndices.end();
    ASSERT_TRUE(isInside);
  }
  for (int idx : atomsOutside) {
    bool isInside = std::find(qmRegionIndices.begin(), qmRegionIndices.end(), idx) != qmRegionIndices.end();
    ASSERT_FALSE(isInside);
  }
}

TEST_F(QmRegionSelectionTests, MultipleQmRegionCandidatesAreConstructedCorrectly) {
  QmRegionSelector qmRegionSelector;
  qmRegionSelector.setLog(silentLogger);
  qmRegionSelector.settings().modifyIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms, {302});
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::initialRadiusForQmRegionSelection, 5.0);
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::cuttingProbability, 0.9);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMinSize, 60);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMaxSize, 75);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionRefMaxSize, 120);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::numAttemptsPerRadius, 30);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::maxNumRefModels, 7);
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, insulin_connectivity_file);

  // Read in structure
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(insulin_xyz_file).first;
  Eigen::RowVector3d centerAtomPosition = structure.getPosition(302);

  // Get bond orders
  auto listsOfNeighbors = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(insulin_connectivity_file);
  auto bondOrders = SwooseUtilities::TopologyUtils::generateBondOrderMatrixFromListsOfNeighbors(listsOfNeighbors);

  // Generate candidates and reference models
  std::vector<QmmmModel> qmmmModelCandidates, qmmmReferenceModels;
  QmRegionCandidateGenerator::generateQmRegionCandidates(qmmmModelCandidates, qmmmReferenceModels, structure,
                                                         bondOrders, qmRegionSelector.settings(), silentLogger);

  // Asserts
  ASSERT_TRUE(qmmmModelCandidates.size() > 10);
  for (const auto& model : qmmmModelCandidates) {
    ASSERT_THAT(model.structure.size(), Eq(model.qmAtomIndices.size()));
    ASSERT_TRUE(model.structure.size() <= 75 && model.structure.size() >= 60);
    ASSERT_THAT(model.molecularCharge, Eq(0));
    ASSERT_THAT(model.spinMultiplicity, Eq(1));
    ASSERT_TRUE(std::find(model.qmAtomIndices.begin(), model.qmAtomIndices.end(), 302) != model.qmAtomIndices.end());
    Eigen::RowVector3d positionZero = model.structure.getPosition(0);
    ASSERT_THAT((positionZero - centerAtomPosition).norm(), DoubleNear(0.0, 1e-6));
  }
  ASSERT_THAT(qmmmReferenceModels.size(), Eq(7));
  for (const auto& model : qmmmReferenceModels) {
    ASSERT_THAT(model.structure.size(), Eq(model.qmAtomIndices.size()));
    ASSERT_TRUE(model.structure.size() <= 120 && model.structure.size() >= 114);
    ASSERT_THAT(model.molecularCharge, Eq(0));
    ASSERT_THAT(model.spinMultiplicity, Eq(1));
    ASSERT_TRUE(std::find(model.qmAtomIndices.begin(), model.qmAtomIndices.end(), 302) != model.qmAtomIndices.end());
    Eigen::RowVector3d positionZero = model.structure.getPosition(0);
    ASSERT_THAT((positionZero - centerAtomPosition).norm(), DoubleNear(0.0, 1e-6));
  }

  // The following two tests assert whether wrong settings result in exceptions
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMinSize, 80);
  qmmmModelCandidates.clear();
  qmmmReferenceModels.clear();
  EXPECT_THROW(QmRegionCandidateGenerator::generateQmRegionCandidates(qmmmModelCandidates, qmmmReferenceModels, structure,
                                                                      bondOrders, qmRegionSelector.settings(), silentLogger),
               std::runtime_error);

  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionRefMaxSize, 70);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMinSize, 60);
  qmmmModelCandidates.clear();
  qmmmReferenceModels.clear();
  EXPECT_THROW(QmRegionCandidateGenerator::generateQmRegionCandidates(qmmmModelCandidates, qmmmReferenceModels, structure,
                                                                      bondOrders, qmRegionSelector.settings(), silentLogger),
               std::runtime_error);
}

TEST_F(QmRegionSelectionTests, ReferenceModelEqualsFullSystemIfPossible) {
  QmRegionSelector qmRegionSelector;
  qmRegionSelector.setLog(silentLogger);
  qmRegionSelector.settings().modifyIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms, {302});
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::initialRadiusForQmRegionSelection, 5.0);
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::cuttingProbability, 0.9);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMinSize, 70);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMaxSize, 100);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionRefMaxSize, 500);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::numAttemptsPerRadius, 1);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::maxNumRefModels, 10);
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, insulin_connectivity_file);

  // Read in structure
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(insulin_xyz_file).first;

  // Get bond orders
  auto listsOfNeighbors = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(insulin_connectivity_file);
  auto bondOrders = SwooseUtilities::TopologyUtils::generateBondOrderMatrixFromListsOfNeighbors(listsOfNeighbors);

  // Generate candidates and reference models
  std::vector<QmmmModel> qmmmModelCandidates, qmmmReferenceModels;
  QmRegionCandidateGenerator::generateQmRegionCandidates(qmmmModelCandidates, qmmmReferenceModels, structure,
                                                         bondOrders, qmRegionSelector.settings(), silentLogger);

  // Asserts
  ASSERT_FALSE(qmmmModelCandidates.empty());
  for (const auto& model : qmmmModelCandidates)
    ASSERT_TRUE(model.structure.size() <= 100 && model.structure.size() >= 70);
  ASSERT_THAT(qmmmReferenceModels.size(), Eq(1));
  const QmmmModel& refModel = qmmmReferenceModels.at(0);
  ASSERT_THAT(refModel.structure.size(), Eq(structure.size()));
  for (int i = 0; i < structure.size(); ++i)
    ASSERT_TRUE(std::find(refModel.qmAtomIndices.begin(), refModel.qmAtomIndices.end(), i) != refModel.qmAtomIndices.end());
}

TEST_F(QmRegionSelectionTests, QmmmModelAnalysisIsPerformedCorrectly) {
  QmRegionSelector qmRegionSelector;
  qmRegionSelector.setLog(silentLogger);
  qmRegionSelector.settings().modifyIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms, {15});
  /*
   * Set large tolerance such that all the candidates will be considered
   * regardless of the error assigned to them.
   */
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::tolerancePercentageError, 1e6);
  // Read in structure
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;

  std::vector<int> qmAtomIndices(structure.size());
  std::iota(qmAtomIndices.begin(), qmAtomIndices.end(), 0); // Fill with 0, 1, ..., structure.size().

  std::vector<QmmmModel> qmmmModelCandidates(5);
  for (int i = 0; i < 5; ++i) {
    QmmmModel model = {structure, qmAtomIndices, 0, 1};
    qmmmModelCandidates.at(i) = model;
  }

  QmmmData qmmmData = {{}, {}, {}, 3};
  for (int j = 0; j < 8; ++j) {
    qmmmData.symmetryScores.push_back(1 + 0.1 * j);
    if (j == 2)
      qmmmData.linkAtomNumbers.push_back(j);
    else
      qmmmData.linkAtomNumbers.push_back(j + 5);
    ForcesCollection f = Eigen::MatrixXd::Random(structure.size(), 3);
    qmmmData.forces.push_back(f);
  }

  qmmmData.forces.at(5) = 0.99 * qmmmData.forces.at(1);
  qmmmData.forces.at(6) = 0.98 * qmmmData.forces.at(1);
  qmmmData.forces.at(7) = 0.97 * qmmmData.forces.at(1);

  QmmmModelAnalyzer analyzerOne(qmRegionSelector.settings(), silentLogger, qmmmData, structure, qmmmModelCandidates);
  int selectedQmRegionIndex = analyzerOne.getIndexOfOptimalModel();
  ASSERT_THAT(selectedQmRegionIndex, Eq(2)); // because tolerance is large and 2 has the smallest num. of link atoms.

  // If tolerance is lower, then 1 should be the optimal choice
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::tolerancePercentageError, 20.0);
  QmmmModelAnalyzer analyzerTwo(qmRegionSelector.settings(), silentLogger, qmmmData, structure, qmmmModelCandidates);
  selectedQmRegionIndex = analyzerTwo.getIndexOfOptimalModel();
  ASSERT_THAT(selectedQmRegionIndex, Eq(1));

  // If forces of candidate 1 are empty, it should not be chosen, but 2 should be chosen again
  ForcesCollection emptyForces;
  qmmmData.forces.at(1) = emptyForces;
  QmmmModelAnalyzer analyzerThree(qmRegionSelector.settings(), silentLogger, qmmmData, structure, qmmmModelCandidates);
  selectedQmRegionIndex = analyzerThree.getIndexOfOptimalModel();
  ASSERT_THAT(selectedQmRegionIndex, Eq(2));
}

TEST_F(QmRegionSelectionTests, QmmmModelAnalysisFailureTests) {
  QmRegionSelector qmRegionSelector;
  qmRegionSelector.setLog(silentLogger);
  qmRegionSelector.settings().modifyIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms, {15});
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::tolerancePercentageError, 50.0);
  // Read in structure
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;

  std::vector<int> qmAtomIndices(structure.size());
  std::iota(qmAtomIndices.begin(), qmAtomIndices.end(), 0); // Fill with 0, 1, ..., structure.size().

  QmmmData qmmmData = {{}, {}, {}, 3};
  std::vector<QmmmModel> qmmmModelCandidates;

  std::string exceptionString = "";
  try {
    QmmmModelAnalyzer analyzer(qmRegionSelector.settings(), silentLogger, qmmmData, structure, qmmmModelCandidates);
  }
  catch (const std::runtime_error& e) {
    exceptionString = e.what();
  }
  ASSERT_STREQ(exceptionString.c_str(), "The list of QM/MM model candidates that was given to the analyzer is empty.");

  qmmmModelCandidates.resize(5);
  for (int i = 0; i < 5; ++i) {
    QmmmModel model = {structure, qmAtomIndices, 0, 1};
    qmmmModelCandidates.at(i) = model;
  }

  // Fill qmmmData
  for (int j = 0; j < 8; ++j) {
    qmmmData.symmetryScores.push_back(1 + 0.1 * j);
    qmmmData.linkAtomNumbers.push_back(j + 5);
    ForcesCollection f1 = Eigen::MatrixXd::Random(structure.size(), 3);
    ForcesCollection f2;
    if (j < 5)
      qmmmData.forces.push_back(f1);
    else
      qmmmData.forces.push_back(f2);
  }

  try {
    QmmmModelAnalyzer analyzer(qmRegionSelector.settings(), silentLogger, qmmmData, structure, qmmmModelCandidates);
  }
  catch (const std::runtime_error& e) {
    exceptionString = e.what();
  }
  ASSERT_STREQ(exceptionString.c_str(),
               "None of the calculations for the reference models were completed successfully.");

  // Re-fill qmmmData
  for (int j = 0; j < 8; ++j) {
    ForcesCollection f1 = Eigen::MatrixXd::Random(structure.size(), 3);
    ForcesCollection f2;
    if (j < 5)
      qmmmData.forces.at(j) = f2;
    else
      qmmmData.forces.at(j) = f1;
  }

  try {
    QmmmModelAnalyzer analyzer(qmRegionSelector.settings(), silentLogger, qmmmData, structure, qmmmModelCandidates);
  }
  catch (const std::runtime_error& e) {
    exceptionString = e.what();
  }
  ASSERT_STREQ(exceptionString.c_str(),
               "None of the calculations for the candidate models were completed successfully.");
}

TEST_F(QmRegionSelectionTests, TwoDisconnectedQmAtomsCanBeSet) {
  QmRegionSelector qmRegionSelector;
  qmRegionSelector.setLog(silentLogger);
  // Two atoms that are far away from each other
  qmRegionSelector.settings().modifyIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms, {302, 87});
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::initialRadiusForQmRegionSelection, 3.0);
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::cuttingProbability, 1.0);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMinSize, 60);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMaxSize, 75);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionRefMaxSize, 120);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::numAttemptsPerRadius, 30);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::maxNumRefModels, 7);
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, insulin_connectivity_file);

  // Read in structure
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(insulin_xyz_file).first;

  qmRegionSelector.generateQmRegion(structure);

  auto qmRegion = qmRegionSelector.getQmRegionStructure();
  auto qmRegionIndices = qmRegionSelector.getQmRegionIndices();
  auto qmRegionInfo = qmRegionSelector.getQmRegionChargeAndMultiplicity();

  ASSERT_THAT(qmRegion.size(), Eq(45));
  ASSERT_THAT(qmRegionIndices.size(), Eq(45));
  ASSERT_THAT(qmRegionInfo.first, Eq(0));
  ASSERT_THAT(qmRegionInfo.second, Eq(1));

  // The same atom twice, check that there are no duplicates
  qmRegionSelector.settings().modifyIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms, {87, 87});
  qmRegionSelector.generateQmRegion(structure);

  auto qmRegionDupl = qmRegionSelector.getQmRegionStructure();
  auto qmRegionIndicesDupl = qmRegionSelector.getQmRegionIndices();
  auto qmRegionInfoDupl = qmRegionSelector.getQmRegionChargeAndMultiplicity();

  ASSERT_THAT(qmRegionDupl.size(), Eq(23));
  ASSERT_THAT(qmRegionIndicesDupl.size(), Eq(23));
  ASSERT_THAT(qmRegionInfoDupl.first, Eq(0));
  ASSERT_THAT(qmRegionInfoDupl.second, Eq(1));

  //  Two atoms that are close
  qmRegionSelector.settings().modifyIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms, {87, 98});
  qmRegionSelector.generateQmRegion(structure);

  auto qmRegionOverlap = qmRegionSelector.getQmRegionStructure();
  auto qmRegionIndicesOverlap = qmRegionSelector.getQmRegionIndices();
  auto qmRegionInfoOverlap = qmRegionSelector.getQmRegionChargeAndMultiplicity();

  ASSERT_THAT(qmRegionOverlap.size(), Eq(37));
  ASSERT_THAT(qmRegionIndicesOverlap.size(), Eq(37));
  ASSERT_THAT(qmRegionInfoOverlap.first, Eq(0));
  ASSERT_THAT(qmRegionInfoOverlap.second, Eq(1));
}

TEST_F(QmRegionSelectionTests, TwoSubsystemsCanBeMergedToOneQmRegion) {
  // Read in structure
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(insulin_xyz_file).first;

  // set subsystem 1 manually
  std::vector<int> firstSubgroup;
  for (int i = 0; i < 10; i++) {
    firstSubgroup.push_back(i);
  }
  // set subsystem 2 manually
  std::vector<int> secondSubgroup;
  for (int j = 5; j < 15; j++) {
    secondSubgroup.push_back(j);
  }
  auto subsystemGenerator = [&structure](std::vector<int> subsystemIndices) -> Utils::AtomCollection {
    Utils::AtomCollection subsystem;
    for (const auto& index : subsystemIndices) {
      try {
        subsystem.push_back(structure.at(index));
      }
      catch (...) {
        throw std::runtime_error("Invalid subsystem index!");
      }
    }
    return subsystem;
  };

  std::vector<Utils::AtomCollection> subRegions(2);
  subRegions.at(0) = subsystemGenerator(firstSubgroup);
  subRegions.at(1) = subsystemGenerator(secondSubgroup);

  std::vector<std::vector<int>> listOfMappings(2);
  listOfMappings.at(0) = firstSubgroup;
  listOfMappings.at(1) = secondSubgroup;

  std::vector<int> indices;

  auto qmRegion = SwooseUtilities::FragmentationHelper::mergeSubsystems(indices, subRegions, listOfMappings);
  ASSERT_THAT(qmRegion.size(), Eq(15));
  for (int i = 0; i < 15; i++) {
    // the QM region should contain atoms 0-14 from the main structure
    ASSERT_THAT(indices.at(i), Eq(i));
    ASSERT_EQ(qmRegion.at(i).getElementType(), structure.at(i).getElementType());
    ASSERT_EQ(qmRegion.at(i).getPosition(), structure.at(i).getPosition());
  }
}

// Run this test only in release builds
#ifdef NDEBUG
TEST_F(QmRegionSelectionTests, QmRegionSelectionWorksInDirectMode) {
  auto& manager = Core::ModuleManager::getInstance();
  if (!manager.moduleLoaded("MockModule")) {
    manager.load("Test/TestUtilities");
  }

  YAML::Node yamlNode = YAML::LoadFile(qm_region_selection_yaml_file);
  // Read in structure
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(insulin_xyz_file).first;
  QmRegionSelector qmRegionSelector;
  auto calculator = std::make_shared<QmmmCalculator>();

  auto mmCalculator = std::make_shared<MolecularMechanics::SfamMolecularMechanicsCalculator>();
  auto mmCalculatorAsCalculator = std::dynamic_pointer_cast<Core::Calculator>(mmCalculator);
  auto qmCalculator = manager.get<Core::Calculator>("MOCK-QM", "MockModule");
  calculator->setUnderlyingCalculators({qmCalculator, mmCalculatorAsCalculator});
  calculator->settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, insulin_connectivity_file);
  calculator->settings().modifyString(Utils::SettingsNames::parameterFilePath, insulin_parameter_file);
  calculator->settings().modifyBool(Utils::SettingsNames::electrostaticEmbedding, false);

  calculator->setStructure(structure);

  Utils::nodeToSettings(qmRegionSelector.settings(), yamlNode, true);
  // Add three more settings
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::yamlSettingsFilePath, qm_region_selection_yaml_file);
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, insulin_connectivity_file);
  qmRegionSelector.settings().modifyString(Utils::SettingsNames::parameterFilePath, insulin_parameter_file);

  qmRegionSelector.setUnderlyingCalculator(calculator);
  qmRegionSelector.setLog(silentLogger);

  // Generate QM region and gather results
  qmRegionSelector.generateQmRegion(structure);

  auto resultIndices = qmRegionSelector.getQmRegionIndices();
  auto resultStructure = qmRegionSelector.getQmRegionStructure();
  auto resultCharge = qmRegionSelector.getQmRegionChargeAndMultiplicity().first;
  auto resultMultiplicity = qmRegionSelector.getQmRegionChargeAndMultiplicity().second;

  /*
   * Since the QM calculator is just a mock calculator, it doesn't make sense to assert
   * which structure exactly was chosen in the end.
   */
  ASSERT_TRUE(resultIndices.size() > 63 && resultIndices.size() < 75);
  ASSERT_TRUE(resultIndices.size() <= static_cast<unsigned long>(resultStructure.size()));
  ASSERT_EQ(resultCharge, 0);
  ASSERT_EQ(resultMultiplicity, 1);
}
#endif

} // namespace Tests
} // namespace Scine
