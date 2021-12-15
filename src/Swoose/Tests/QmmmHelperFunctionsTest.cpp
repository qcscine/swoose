/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include <Swoose/MolecularMechanics/SFAM/SfamCalculatorSettings.h>
#include <Swoose/MolecularMechanics/SFAM/SfamMolecularMechanicsCalculator.h>
#include <Swoose/QMMM/QmmmCalculatorSettings.h>
#include <Swoose/QMMM/QmmmGradientsEvaluator.h>
#include <Swoose/QMMM/QmmmHelpers.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <gmock/gmock.h>
#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <numeric>

using namespace testing;
namespace Scine {
using namespace Qmmm;
namespace Tests {

/**
 * @class QmmmHelperFunctionsTest QmmmHelperFunctionsTest.cpp
 * @brief Tests the functionality of some helper functions which are applied in the QM/MM calculator.
 * @test
 */
class QmmmHelperFunctionsTest : public Test {
 public:
  MolecularMechanics::SfamMolecularMechanicsCalculator calculator;
  static constexpr const char* qmRegionXyzFile = "created_qm_region_for_test.xyz";
  static constexpr const char* pointChargesFile = "point_charges_test.pc";
  Utils::AtomCollection testStructure;
  Core::Log silentLogger = Core::Log::silent();

  void SetUp() override {
    testStructure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;
    calculator.setLog(silentLogger);
  }

  ~QmmmHelperFunctionsTest() override {
    boost::filesystem::remove(qmRegionXyzFile);
    boost::filesystem::remove(pointChargesFile);
  }
};

TEST_F(QmmmHelperFunctionsTest, QmRegionIsCorrectlyCreated) {
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  calculator.setStructure(testStructure);

  auto listsOfNeighbors = calculator.listsOfNeighbors();
  std::vector<int> qmAtoms = {30, 31, 32, 23, 22, 29, 21, 28, 20, 26, 27, 19, 24, 25};
  std::vector<int> mmBoundaryAtoms;

  auto qmRegion = QmmmHelpers::createQmRegion(qmAtoms, testStructure, listsOfNeighbors, qmRegionXyzFile, mmBoundaryAtoms);

  ASSERT_THAT(qmRegion.size() - 1, Eq(qmAtoms.size())); // -1 because of the link atom
  ASSERT_THAT(mmBoundaryAtoms.size(), Eq(1));

  Eigen::RowVector3d pos1 = testStructure.getPosition(qmAtoms.at(5));
  Eigen::RowVector3d pos2 = qmRegion.getPosition(5);
  Eigen::RowVector3d pos3 = testStructure.getPosition(qmAtoms.at(11));
  Eigen::RowVector3d pos4 = qmRegion.getPosition(11);
  for (int k = 0; k < 3; ++k) {
    ASSERT_THAT(pos1(k), DoubleNear(pos2(k), 1e-6));
    ASSERT_THAT(pos3(k), DoubleNear(pos4(k), 1e-6));
  }

  // Read in XYZ file of QM region
  Utils::AtomCollection qmRegionFromFile = Utils::ChemicalFileHandler::read(qmRegionXyzFile).first;
  ASSERT_THAT(qmRegion.size(), Eq(qmRegionFromFile.size()));

  for (int i = 0; i < qmRegion.size(); ++i) {
    ASSERT_THAT(qmRegion.getElement(i), Eq(qmRegionFromFile.getElement(i)));
    ASSERT_THAT(qmRegion.getPosition(i).x(), DoubleNear(qmRegionFromFile.getPosition(i).x(), 1e-6));
    ASSERT_THAT(qmRegion.getPosition(i).y(), DoubleNear(qmRegionFromFile.getPosition(i).y(), 1e-6));
    ASSERT_THAT(qmRegion.getPosition(i).z(), DoubleNear(qmRegionFromFile.getPosition(i).z(), 1e-6));
  }
}

TEST_F(QmmmHelperFunctionsTest, LinkAtomIsCorrectlyAdded) {
  Eigen::RowVector3d qmAtomPos(0.0, 0.0, 0.0);
  Utils::Atom qmAtom(Utils::ElementType::C, qmAtomPos);

  Eigen::RowVector3d mmAtomPos(1.5, 0.0, 0.0);
  Utils::Atom mmAtom(Utils::ElementType::O, mmAtomPos);

  Utils::AtomCollection qmRegion;
  qmRegion.push_back(qmAtom);

  QmmmHelpers::addOneLinkAtom(qmRegion, qmAtom, mmAtom);

  ASSERT_THAT(qmRegion.size(), Eq(2));
  ASSERT_THAT(qmRegion.getPosition(0).x(), DoubleNear(qmAtom.getPosition().x(), 1e-6));

  Eigen::RowVector3d expectedLinkAtomPos(Utils::ElementInfo::covalentRadius(Utils::ElementType::C) +
                                             Utils::ElementInfo::covalentRadius(Utils::ElementType::H),
                                         0.0, 0.0);

  ASSERT_THAT(qmRegion.getElement(1), Eq(Utils::ElementType::H));
  for (int k = 0; k < 3; ++k) {
    ASSERT_THAT(expectedLinkAtomPos(k), DoubleNear(qmRegion.getPosition(1)(k), 1e-6));
  }
}

TEST_F(QmmmHelperFunctionsTest, PointChargesFileIsCorrectlyWritten) {
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  calculator.setStructure(testStructure);

  auto listsOfNeighbors = calculator.listsOfNeighbors();
  std::vector<int> qmAtoms = {30, 31, 32, 23, 22, 29, 21, 28, 20, 26, 27, 19, 24, 25};

  QmmmHelpers::ChargeRedistributionResult redistributedChargesResult;
  redistributedChargesResult.atomicCharges = calculator.atomicCharges();
  redistributedChargesResult.auxiliaryCharges.push_back(1.0);
  redistributedChargesResult.auxiliaryCharges.push_back(2.0);
  redistributedChargesResult.auxiliaryCharges.push_back(3.0);
  redistributedChargesResult.positionsOfAuxiliaryCharges.resize(3, 3);
  auto pos1 = Eigen::RowVector3d(0.0, 0.0, 5.0);
  auto pos2 = Eigen::RowVector3d(0.0, 2.0, 3.0);
  auto pos3 = Eigen::RowVector3d(1.0, 2.0, 3.0);
  redistributedChargesResult.positionsOfAuxiliaryCharges.row(0) = pos1;
  redistributedChargesResult.positionsOfAuxiliaryCharges.row(1) = pos2;
  redistributedChargesResult.positionsOfAuxiliaryCharges.row(2) = pos3;

  QmmmHelpers::writePointChargesFile(calculator.getPositions(), redistributedChargesResult, qmAtoms, pointChargesFile);

  ASSERT_TRUE(boost::filesystem::exists(pointChargesFile));

  std::string content1;
  std::string content2;
  std::string content3;
  std::ifstream check;
  check.open(pointChargesFile);
  if (check.is_open()) {
    check >> content1;
    check >> content2;
    check >> content3;
  }
  check.close();

  std::string expectedContent =
      std::to_string(testStructure.size() - qmAtoms.size() + redistributedChargesResult.auxiliaryCharges.size());
  ASSERT_THAT(content1, Eq(expectedContent));
  ASSERT_THAT(redistributedChargesResult.auxiliaryCharges.size(),
              Eq(redistributedChargesResult.positionsOfAuxiliaryCharges.rows()));

  auto content2AsDouble = std::stod(content2);
  auto content3AsDouble = std::stod(content3);
  ASSERT_THAT(content2AsDouble, DoubleNear(calculator.atomicCharges().at(0), 1e-6));
  ASSERT_THAT(content3AsDouble * Utils::Constants::bohr_per_angstrom, DoubleNear(testStructure.getPosition(0).x(), 1e-6));
}

TEST_F(QmmmHelperFunctionsTest, QmmmRegionValidityCheckWorks) {
  std::vector<int> qmAtoms = {30, 31, 32, 23, 22, 29, 21, 28, 20, 26, 27, 19, 24, 25}; // this vector should be valid
  EXPECT_NO_THROW(QmmmHelpers::checkValidityOfQmRegion(qmAtoms, testStructure));

  qmAtoms = {30, 31, 32, 23, 22, 29, 21, -1, 20, 26, 27, 19, 24, 25}; // index -1 is not valid
  EXPECT_THROW(QmmmHelpers::checkValidityOfQmRegion(qmAtoms, testStructure), std::runtime_error);

  qmAtoms = {30, 31, 32, 23, 1000, 29, 21, 28, 20, 26, 27, 19, 24, 25}; // index 1000 is not valid
  EXPECT_THROW(QmmmHelpers::checkValidityOfQmRegion(qmAtoms, testStructure), std::runtime_error);

  qmAtoms = {30, 31, 32, 23, 33, 29, 21, 28, 20, 26, 27, 19, 24, 25}; // index 33 is not valid
  EXPECT_THROW(QmmmHelpers::checkValidityOfQmRegion(qmAtoms, testStructure), std::runtime_error);

  qmAtoms = {30, 31, 32, 23, 22, 29, 21, 28, 20, 26, 27, 31, 24, 25}; // Index 31 exists twice
  EXPECT_THROW(QmmmHelpers::checkValidityOfQmRegion(qmAtoms, testStructure), std::runtime_error);
}

TEST_F(QmmmHelperFunctionsTest, TotalGradientsAreCorrectlyAssembled) {
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  // Basically no cutoff radius:
  calculator.settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 20.0);

  calculator.setStructure(testStructure);
  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);

  auto results = calculator.calculate("test calculation 1");
  const Utils::GradientCollection& mmGradients = results.get<Utils::Property::Gradients>();

  std::vector<int> qmAtoms = {30, 31, 32, 23};
  Utils::GradientCollection qmGradients(5, 3); // last atom shall be link atom
  qmGradients << 100, 100, 100, 200, 200, 200, 300, 300, 300, 400, 400, 400, 500, 500, 500;

  std::vector<int> mmBoundaryAtoms;
  auto qmRegion = QmmmHelpers::createQmRegion(qmAtoms, testStructure, calculator.listsOfNeighbors(), "", mmBoundaryAtoms);

  ASSERT_THAT(mmBoundaryAtoms.size(), Eq(1));
  ASSERT_THAT(mmBoundaryAtoms.at(0), Eq(22));

  // Artificial gradients contributions from the point charges.
  Utils::GradientCollection pcGradients(30, 3);
  for (int k = 0; k < 30; ++k) {
    for (int l = 0; l < 3; ++l) {
      pcGradients(k, l) = 0.0;
    }
  }

  QmmmGradientsEvaluator gradEvaluator(qmGradients, mmGradients, pcGradients, qmAtoms, mmBoundaryAtoms,
                                       calculator.listsOfNeighbors(), testStructure, qmRegion);
  Utils::GradientCollection combinedGradients = gradEvaluator.calculateQmmmGradients();

  ASSERT_THAT(combinedGradients.rows(), Eq(mmGradients.rows()));

  for (int i = 0; i < combinedGradients.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i == 30) {
        ASSERT_THAT(combinedGradients(i, j), DoubleNear(qmGradients(0, j), 1.0));
      }
      else if (i == 31) {
        ASSERT_THAT(combinedGradients(i, j), DoubleNear(qmGradients(1, j), 1.0));
      }
      else if (i == 32) {
        ASSERT_THAT(combinedGradients(i, j), DoubleNear(qmGradients(2, j), 1.0));
      }
      else if (i == 23) {
        ASSERT_FALSE(std::abs(combinedGradients(i, j) - qmGradients(3, j)) < 1.0);
      }
      else if (i == 22) {
        ASSERT_FALSE(std::abs(combinedGradients(i, j) - mmGradients(i, j)) < 1.0);
      }

      else {
        ASSERT_THAT(combinedGradients(i, j), DoubleNear(mmGradients(i, j), 1e-6));
      }
    }
  }
}

TEST_F(QmmmHelperFunctionsTest, ChargesAreCorrectlyRedistributed) {
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  calculator.setStructure(testStructure);

  std::vector<int> qmAtoms = {30, 31, 32, 23, 22, 29, 21, 28, 20, 26, 27, 19, 24, 25};
  std::vector<int> mmBoundaryAtoms = {2};

  // First test the simple redistribution scheme with no auxiliary point charges:
  std::string redistributedChargeScheme = SwooseUtilities::OptionNames::redistributedChargeOption;
  auto redistributedChargesResult =
      QmmmHelpers::getRedistributedCharges(calculator.atomicCharges(), testStructure.getPositions(), mmBoundaryAtoms,
                                           calculator.listsOfNeighbors(), qmAtoms, redistributedChargeScheme);
  auto& newCharges = redistributedChargesResult.atomicCharges;

  ASSERT_THAT(newCharges.size(), Eq(calculator.atomicCharges().size()));

  for (int i = 0; i < newCharges.size(); ++i) {
    if (i == 2) {
      ASSERT_THAT(newCharges.at(i), DoubleNear(0.0, 1e-6));
    }
    else if (i == 1) {
      double expectedCharge = calculator.atomicCharges().at(i) + calculator.atomicCharges().at(2) / 2;
      ASSERT_THAT(newCharges.at(i), DoubleNear(expectedCharge, 1e-6));
    }
    else if (i == 3) {
      double expectedCharge = calculator.atomicCharges().at(i) + calculator.atomicCharges().at(2) / 2;
      ASSERT_THAT(newCharges.at(i), DoubleNear(expectedCharge, 1e-6));
    }
    else {
      ASSERT_THAT(newCharges.at(i), DoubleNear(calculator.atomicCharges().at(i), 1e-6));
    }
    // Make sure charges were not all zero due to some bug:
    ASSERT_TRUE(std::abs(calculator.atomicCharges().at(i)) > 1e-2);
  }

  // Check that the sum of charges is still the same
  double expectedSum = std::accumulate(calculator.atomicCharges().begin(), calculator.atomicCharges().end(), 0.0);
  double sum = std::accumulate(newCharges.begin(), newCharges.end(), 0.0);
  ASSERT_THAT(sum, DoubleNear(expectedSum, 1e-6));

  // Now test the charge- and dipole-conserving redistribution scheme with auxiliary point charges on bond vectors:
  redistributedChargeScheme = SwooseUtilities::OptionNames::redistributedChargeAndDipolesOption;
  redistributedChargesResult =
      QmmmHelpers::getRedistributedCharges(calculator.atomicCharges(), testStructure.getPositions(), mmBoundaryAtoms,
                                           calculator.listsOfNeighbors(), qmAtoms, redistributedChargeScheme);
  auto& newChargesDipoleScheme = redistributedChargesResult.atomicCharges;
  auto& auxCharges = redistributedChargesResult.auxiliaryCharges;
  auto& auxPositions = redistributedChargesResult.positionsOfAuxiliaryCharges;

  ASSERT_THAT(newChargesDipoleScheme.size(), Eq(calculator.atomicCharges().size()));
  ASSERT_THAT(auxCharges.size(), Eq(2));
  ASSERT_THAT(auxCharges.size(), Eq(auxPositions.rows()));

  for (int i = 0; i < newChargesDipoleScheme.size(); ++i) {
    if (i == 2) {
      ASSERT_THAT(newChargesDipoleScheme.at(i), DoubleNear(0.0, 1e-6));
    }
    else if (i == 1) {
      double expectedCharge = calculator.atomicCharges().at(i) - calculator.atomicCharges().at(2) / 2;
      ASSERT_THAT(newChargesDipoleScheme.at(i), DoubleNear(expectedCharge, 1e-6));
    }
    else if (i == 3) {
      double expectedCharge = calculator.atomicCharges().at(i) - calculator.atomicCharges().at(2) / 2;
      ASSERT_THAT(newChargesDipoleScheme.at(i), DoubleNear(expectedCharge, 1e-6));
    }
    else {
      ASSERT_THAT(newChargesDipoleScheme.at(i), DoubleNear(calculator.atomicCharges().at(i), 1e-6));
    }
  }

  // Check the auxiliary virtual charges
  ASSERT_THAT(auxCharges.at(0), DoubleNear(calculator.atomicCharges().at(2), 1e-6));
  ASSERT_THAT(auxCharges.at(1), DoubleNear(calculator.atomicCharges().at(2), 1e-6));
  Eigen::RowVector3d posOfAuxCharge0 = 0.5 * (testStructure.getPosition(2) + testStructure.getPosition(1));
  Eigen::RowVector3d posOfAuxCharge1 = 0.5 * (testStructure.getPosition(2) + testStructure.getPosition(3));

  for (int k = 0; k < 3; ++k) {
    ASSERT_THAT(auxPositions(0, k), DoubleNear(posOfAuxCharge0(k), 1e-6));
    ASSERT_THAT(auxPositions(1, k), DoubleNear(posOfAuxCharge1(k), 1e-6));
  }

  // Check that the sum of charges is still the same
  double newSum = std::accumulate(newChargesDipoleScheme.begin(), newChargesDipoleScheme.end(), 0.0);
  for (const auto& aux : auxCharges)
    newSum += aux;
  ASSERT_THAT(newSum, DoubleNear(expectedSum, 1e-6));
}

} // namespace Tests
} // namespace Scine
