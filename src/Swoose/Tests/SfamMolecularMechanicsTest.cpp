/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include <Swoose/MolecularMechanics/SFAM/SfamCalculatorSettings.h>
#include <Swoose/MolecularMechanics/SFAM/SfamMolecularMechanicsCalculator.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <gmock/gmock.h>
#include <boost/filesystem/operations.hpp>

using namespace testing;
namespace Scine {
using namespace MolecularMechanics;
using namespace SwooseUtilities;
namespace Tests {

/**
 * @class ASfamMolecularMechanicsTest SfamMolecularMechanicsTest.cpp
 * @brief Tests concerning SFAM's molecular mechanics method in SCINE.
 * @test
 */
class ASfamMolecularMechanicsTest : public Test {
 public:
  SfamMolecularMechanicsCalculator calculator;
  Core::Log silentLogger = Core::Log::silent();
  void SetUp() override {
    calculator.setLog(silentLogger);
  }
};

TEST_F(ASfamMolecularMechanicsTest, SettingsAreSetCorrectly) {
  using namespace SwooseUtilities::SettingsNames;
  calculator.settings().modifyBool(onlyCalculateBondedContribution, true);
  calculator.settings().modifyBool(printContributionsMolecularMechanics, true);
  calculator.settings().modifyString(Utils::SettingsNames::parameterFilePath, "pfp");
  calculator.settings().modifyString(connectivityFilePath, "cfp");
  calculator.settings().modifyString(sfamAtomTypeLevel, "unique");
  calculator.settings().modifyBool(detectBondsWithCovalentRadii, false);
  calculator.settings().modifyDouble(nonCovalentCutoffRadius, 4.2);
  calculator.settings().modifyBool(hydrogenBondCorrection, false);

  ASSERT_THAT(calculator.settings().getBool(onlyCalculateBondedContribution), Eq(true));
  ASSERT_THAT(calculator.settings().getBool(printContributionsMolecularMechanics), Eq(true));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::parameterFilePath), Eq("pfp"));
  ASSERT_THAT(calculator.settings().getString(connectivityFilePath), Eq("cfp"));
  ASSERT_THAT(calculator.settings().getString(sfamAtomTypeLevel), Eq("unique"));
  ASSERT_THAT(calculator.settings().getBool(detectBondsWithCovalentRadii), Eq(false));
  ASSERT_THAT(calculator.settings().getDouble(nonCovalentCutoffRadius), Eq(4.2));
  ASSERT_THAT(calculator.settings().getBool(hydrogenBondCorrection), Eq(false));
}

// The energy of the test structure that is calculated here contains contributions of all force field interaction types.
TEST_F(ASfamMolecularMechanicsTest, EnergyAndGradientsAreCalculatedCorrectlyForMelatonin) {
  Utils::AtomCollection testStructure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;

  calculator.settings().modifyString(Utils::SettingsNames::parameterFilePath, melatonin_param_file);
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  calculator.settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 20.0); // basically no
                                                                                                     // cutoff radius

  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);

  calculator.setStructure(testStructure);

  // Clone calculator to test clone interface as well
  auto clonedCalculator = calculator.clone();

  // Calculate and gather results
  Utils::Results results = clonedCalculator->calculate("test calculation 1");
  double energy = results.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;
  Utils::GradientCollection gradients = results.get<Utils::Property::Gradients>() * Utils::Constants::kCalPerMol_per_hartree;
  Utils::HessianMatrix hessian = results.get<Utils::Property::Hessian>() * Utils::Constants::kCalPerMol_per_hartree;

  // Numerical gradients
  Utils::PositionCollection positions = testStructure.getPositions();
  Utils::GradientCollection numGrad(positions.rows(), 3);

  double stepsize = 1e-5;

  for (int i = 0; i < positions.rows(); ++i) {
    Eigen::RowVector3d pos = positions.row(i);
    Eigen::RowVector3d tmpPos = pos;
    for (int j = 0; j < 3; ++j) {
      pos[j] += stepsize;
      positions.row(i) = pos;

      clonedCalculator->setRequiredProperties(Utils::Property::Energy); // This is needed, but why? Bug?
      clonedCalculator->modifyPositions(positions);
      results = clonedCalculator->calculate("");
      double energyPlus = results.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;

      pos = tmpPos;
      pos[j] -= stepsize;
      positions.row(i) = pos;

      clonedCalculator->modifyPositions(positions);
      results = clonedCalculator->calculate("");
      double energyMinus = results.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;

      pos = tmpPos;
      positions.row(i) = pos;

      numGrad(i, j) = (energyPlus - energyMinus) / (2.0 * stepsize);
    }
  }

  // Assertions
  ASSERT_THAT(energy, DoubleNear(-4.009964844, 1e-6));
  ASSERT_THAT(hessian(0, 0), DoubleNear(92.95714161, 1e-3));
  ASSERT_THAT(hessian(2, 7), DoubleNear(0.6978703803, 1e-3));
  ASSERT_THAT(hessian(11, 11), DoubleNear(381.9659289, 1e-3));

  for (int k = 0; k < gradients.rows(); ++k) {
    for (int l = 0; l < 3; ++l) {
      ASSERT_THAT(gradients(k, l), DoubleNear(numGrad(k, l), 1e-3));
    }
  }

  // Check if SFAM is translationally invariant
  Utils::PositionCollection translatedPositions(positions.rows(), 3);
  for (int i = 0; i < positions.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      translatedPositions.row(i)[j] = positions.row(i)[j] + 1.0;
    }
  }
  calculator.modifyPositions(translatedPositions);
  Utils::Results resultsAfterTranslation = calculator.calculate("for translated molecule");
  double transEnergy = resultsAfterTranslation.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;
  auto transGradients = resultsAfterTranslation.get<Utils::Property::Gradients>() * Utils::Constants::kCalPerMol_per_hartree;
  ASSERT_THAT(transEnergy, DoubleNear(energy, 1e-6));
  for (int k = 0; k < transGradients.rows(); ++k) {
    for (int l = 0; l < 3; ++l) {
      ASSERT_THAT(transGradients(k, l), DoubleNear(gradients(k, l), 1e-3));
    }
  }

  calculator.modifyPositions(positions);
  /*
   * The same calculation with the original calculator and no connectivity file given.
   */
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, "");
  Utils::Results newResults = calculator.calculate("test calculation 2");
  double newEnergy = newResults.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;
  auto newGradients = newResults.get<Utils::Property::Gradients>() * Utils::Constants::kCalPerMol_per_hartree;
  auto newHessian = newResults.get<Utils::Property::Hessian>() * Utils::Constants::kCalPerMol_per_hartree;
  auto success = results.get<Utils::Property::SuccessfulCalculation>();

  // Assertions
  ASSERT_THAT(success, Eq(true));
  ASSERT_THAT(energy, DoubleNear(newEnergy, 1e-6));
  ASSERT_THAT(hessian(0, 0), DoubleNear(newHessian(0, 0), 1e-3));
  ASSERT_THAT(hessian(2, 7), DoubleNear(newHessian(2, 7), 1e-3));
  ASSERT_THAT(hessian(11, 11), DoubleNear(newHessian(11, 11), 1e-3));
  for (int k = 0; k < gradients.rows(); ++k) {
    for (int l = 0; l < 3; ++l) {
      ASSERT_THAT(gradients(k, l), DoubleNear(newGradients(k, l), 1e-3));
    }
  }

  /*
   * The same calculation with a parameter file that is missing C6 coefficients.
   * The energy should be similar, but not quite the same as the C6 now deviate slightly
   * from the other calculation.
   */
  constexpr const char* tmpFile = "tmp_parameter_file.dat";
  auto silentLogger = Core::Log::silent();
  calculator.setLog(silentLogger);
  Utils::FilesystemHelpers::copyFile(melatonin_param_file_no_c6, tmpFile);
  calculator.settings().modifyString(Utils::SettingsNames::parameterFilePath, tmpFile);
  calculator.setRequiredProperties(Utils::Property::Energy);
  calculator.setStructure(testStructure); // new set structure because re-initialization is necessary
  newResults = calculator.calculate("test calculation 3");

  // Clean-up
  boost::filesystem::remove_all(tmpFile);

  // Assertions
  double otherC6Energy = newResults.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;
  ASSERT_TRUE(std::abs(energy - otherC6Energy) > 1e-3);
  ASSERT_TRUE(std::abs(energy - otherC6Energy) < 1e-1);
}

// Tests that the Hessian is identical to a fully numerical one (Hessian is half numerical and half analytical).
TEST_F(ASfamMolecularMechanicsTest, HessianIsIdenticalToTheFullyNumericalHessianForMelatonin) {
  Utils::AtomCollection testStructure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;

  calculator.settings().modifyString(Utils::SettingsNames::parameterFilePath, melatonin_param_file);
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  calculator.settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 20.0); // basically no
                                                                                                     // cutoff radius

  calculator.setRequiredProperties(Utils::Property::Hessian);
  calculator.setStructure(testStructure);

  // Calculate and gather results
  Utils::Results results = calculator.calculate("test calculation 3");
  const Utils::HessianMatrix& hessian = results.get<Utils::Property::Hessian>();

  // Calculated Hessian fully numerically
  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  Utils::NumericalHessianCalculator numericalHessianCalculator(calculator);
  Utils::Results numericalResult = numericalHessianCalculator.calculate();
  const Utils::HessianMatrix& numericalHessian = numericalResult.get<Utils::Property::Hessian>();

  // Assertions
  ASSERT_THAT(hessian.rows(), Eq(numericalHessian.rows()));
  ASSERT_THAT(hessian.cols(), Eq(numericalHessian.cols()));
  for (int i = 0; i < hessian.rows(); ++i) {
    for (int j = 0; j < hessian.cols(); ++j) {
      ASSERT_THAT(hessian(i, j), DoubleNear(numericalHessian(i, j), 1e-6));
    }
  }
}

TEST_F(ASfamMolecularMechanicsTest, TopologyUtilsDividesAStructureCorrectly) {
  // Does it correctly divide two bonded atoms?
  auto bo1 = Utils::BondOrderCollection(2);
  bo1.setZero();
  bo1.setOrder(1, 0, 1.0);

  bool ok1;
  auto r1 = TopologyUtils::divideStructureAtBond(0, 1, bo1, ok1);
  ASSERT_TRUE(ok1);
  ASSERT_THAT(r1[0], Eq(1));
  ASSERT_THAT(r1[1], Eq(2));

  // Does it correctly divide three bonded atoms?
  auto bo2 = Utils::BondOrderCollection(3);
  bo2.setZero();
  bo2.setOrder(1, 0, 1.0);
  bo2.setOrder(2, 1, 1.0);

  bool ok2;
  auto r2 = TopologyUtils::divideStructureAtBond(0, 1, bo2, ok2);
  ASSERT_TRUE(ok2);
  ASSERT_THAT(r2[0], Eq(1));
  ASSERT_THAT(r2[1], Eq(2));
  ASSERT_THAT(r2[2], Eq(2));

  // Does it detect a cycle?
  auto bo3 = Utils::BondOrderCollection(3);
  bo3.setZero();
  bo3.setOrder(1, 0, 1.0);
  bo3.setOrder(2, 0, 1.0);
  bo3.setOrder(2, 1, 1.0);

  bool ok3;
  auto r3 = TopologyUtils::divideStructureAtBond(0, 1, bo3, ok3);
  ASSERT_FALSE(ok3);
}

TEST_F(ASfamMolecularMechanicsTest, IncorrectConnectivityFileCausesException) {
  Utils::AtomCollection testStructure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;
  calculator.settings().modifyString(Utils::SettingsNames::parameterFilePath, melatonin_param_file);
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath,
                                     melatonin_connectivity_file_with_error);
  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);

  std::string exceptionString = "";
  try {
    calculator.setStructure(testStructure);
  }
  catch (const std::runtime_error& e) {
    exceptionString = e.what();
  }

  ASSERT_STREQ(exceptionString.c_str(), "Connectivity file is invalid! Error during check of atom with index: 32");
}

} // namespace Tests
} // namespace Scine
