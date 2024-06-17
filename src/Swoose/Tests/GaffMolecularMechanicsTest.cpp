/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include <Swoose/MolecularMechanics/GAFF/GaffAtomTypeIdentifier.h>
#include <Swoose/MolecularMechanics/GAFF/GaffCalculatorSettings.h>
#include <Swoose/MolecularMechanics/GAFF/GaffMolecularMechanicsCalculator.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondDetector.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>
#include <fstream>

using namespace testing;
namespace Scine {
using namespace MolecularMechanics;
namespace Tests {

/**
 * @class AGaffMolecularMechanicsTest GaffMolecularMechanicsTest.cpp
 * @brief Tests concerning the GAFF molecular mechanics method in SCINE.
 * @test
 */
class AGaffMolecularMechanicsTest : public Test {
 public:
  GaffMolecularMechanicsCalculator calculator;
  std::string atomicChargesFile = "atomic_charges_for_gaff_test.dat";
  Core::Log silentLogger = Core::Log::silent();
  void SetUp() override {
    calculator.setLog(silentLogger);
  }
  void TearDown() override {
    boost::filesystem::remove_all(atomicChargesFile);
  }
};

TEST_F(AGaffMolecularMechanicsTest, GaffAtomTypesAreAssignedCorrectlyForAlanine) {
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;
  Utils::BondOrderCollection bondOrders = Utils::BondDetector::detectBonds(structure);
  auto listsOfNeighbors =
      SwooseUtilities::TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(structure.size(), bondOrders, 0.5);

  GaffAtomTypeIdentifier gaffAtomTypeIdentifier(structure.size(), structure.getElements(), listsOfNeighbors);
  auto atomTypes = gaffAtomTypeIdentifier.getAtomTypes();

  const char* const correctTypes[] = {"hn", "n3", "hn", "c3", "hc", "c3", "hc", "hc", "hc", "c", "o", "oh", "ho"};
  for (int i = 0; i < atomTypes.size(); ++i) {
    ASSERT_STREQ(atomTypes.getAtomType(i).c_str(), correctTypes[i]);
  }
}

TEST_F(AGaffMolecularMechanicsTest, GaffAtomTypesAreAssignedCorrectlyForMelatonin) {
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;
  Utils::BondOrderCollection bondOrders = Utils::BondDetector::detectBonds(structure);
  auto listsOfNeighbors =
      SwooseUtilities::TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(structure.size(), bondOrders, 0.5);

  GaffAtomTypeIdentifier gaffAtomTypeIdentifier(structure.size(), structure.getElements(), listsOfNeighbors, "");
  auto atomTypes = gaffAtomTypeIdentifier.getAtomTypes();

  const char* const correctTypes[] = {"?",  "?",  "?",  "ca", "ca", "ca", "ca", "ca", "ca", "hn", "ha",
                                      "ha", "ha", "ha", "os", "c3", "hc", "hc", "hc", "c3", "c3", "n",
                                      "c",  "c3", "hc", "hc", "hc", "hc", "hn", "o",  "hc", "hc", "hc"};
  // TODO: Atoms 0, 1 and 2 are difficult to assign.
  for (int i = 3; i < atomTypes.size(); ++i) {
    ASSERT_STREQ(atomTypes.getAtomType(i).c_str(), correctTypes[i]);
  }
}

TEST_F(AGaffMolecularMechanicsTest, GaffDummyAtomTypesAreReadFromFileCorrectlyForMelatonin) {
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;
  Utils::BondOrderCollection bondOrders = Utils::BondDetector::detectBonds(structure);
  auto listsOfNeighbors =
      SwooseUtilities::TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(structure.size(), bondOrders, 0.5);

  GaffAtomTypeIdentifier gaffAtomTypeIdentifier(structure.size(), structure.getElements(), listsOfNeighbors,
                                                gaff_dummy_atom_types_file);
  auto atomTypes = gaffAtomTypeIdentifier.getAtomTypes();

  const char* const correctTypes[] = {"x1", "x2", "x3", "x4", "x5", "ca", "ca", "ca", "ca", "hn", "ha",
                                      "ha", "ha", "ha", "os", "c3", "xx", "yy", "zz", "c3", "c3", "n",
                                      "c",  "c3", "hc", "hc", "hc", "hc", "y5", "y4", "y3", "y2", "y1"};
  // Some atoms are replaced by non-existing dummy atom types (containing x, y, or z), however, should be
  // fine, because the atom types are read in from file.
  for (int i = 0; i < atomTypes.size(); ++i) {
    ASSERT_STREQ(atomTypes.getAtomType(i).c_str(), correctTypes[i]);
  }
}

TEST_F(AGaffMolecularMechanicsTest, GaffAtomTypesAreAssignedCorrectlyForTyrosine) {
  std::stringstream ss("24\n\n"
                       "N      1.1085     -2.1880      0.3267\n"
                       "C      1.4287     -0.9918     -0.4835\n"
                       "C      2.8837     -0.9704     -0.9469\n"
                       "O      3.8526     -1.5054     -0.4376\n"
                       "O      3.1229     -0.2719     -2.0788\n"
                       "C      1.1432      0.3020      0.2994\n"
                       "H      1.6116      0.2640      1.3046\n"
                       "H      1.6185      1.1584     -0.2213\n"
                       "C     -0.3239      0.5597      0.4153\n"
                       "C     -1.0257      1.0790     -0.6765\n"
                       "H     -0.4944      1.2871     -1.6127\n"
                       "C     -2.3860      1.3333     -0.5880\n"
                       "H     -2.9335      1.7409     -1.4450\n"
                       "C     -3.0526      1.0606      0.6153\n"
                       "C     -2.3612      0.5402      1.7163\n"
                       "H     -2.8805      0.3238      2.6568\n"
                       "C     -0.9982      0.2932      1.6064\n"
                       "H     -0.4531     -0.1198      2.4631\n"
                       "O     -4.3928      1.3349      0.6371\n"
                       "H     -4.7252      1.0949      1.4930\n"
                       "H      0.7619     -1.0294     -1.3829\n"
                       "H      4.0509     -0.3011     -2.2915\n"
                       "H      1.6183     -2.1756      1.1852\n"
                       "H      1.3339     -3.0135     -0.1856\n");

  Utils::AtomCollection structure = Utils::XyzStreamHandler::read(ss);
  Utils::BondOrderCollection bondOrders = Utils::BondDetector::detectBonds(structure);
  auto listsOfNeighbors =
      SwooseUtilities::TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(structure.size(), bondOrders, 0.5);

  // "correct" mean that these are the ones deducted by OpenBabel
  const char* const correctTypes[] = {"n3", "c3", "c",  "o",  "oh", "c3", "hc", "hc", "ca", "ca", "ha", "ca",
                                      "ha", "ca", "ca", "ha", "ca", "ha", "oh", "ho", "h1", "ho", "hn", "hn"};
  GaffAtomTypeIdentifier gaffAtomTypeIdentifier(structure.size(), structure.getElements(), listsOfNeighbors);
  auto atomTypes = gaffAtomTypeIdentifier.getAtomTypes();
  for (int i = 0; i < atomTypes.size(); ++i) {
    if (strcmp(correctTypes[i], "h1") == 0) {
      // "hc" is the "non-special" equivalent atom type
      ASSERT_TRUE(atomTypes.getAtomType(i) == std::string("h1") || atomTypes.getAtomType(i) == std::string("hc"));
      continue;
    }
    ASSERT_STREQ(atomTypes.getAtomType(i).c_str(), correctTypes[i]);
  }
}

TEST_F(AGaffMolecularMechanicsTest, GaffCalculationOfDistortedEthaneMolecule) {
  std::stringstream ss("8\n\n"
                       "C 0.683716 4.459351 0.439667\n"
                       "H 0.786743 4.051086 1.451935\n"
                       "H 0.160339 3.709412 -0.165070\n"
                       "H -0.056887 5.433069 0.544398\n"
                       "C 2.036743 4.810486 -0.158880\n"
                       "H 2.823608 4.143188 0.212402\n"
                       "H 2.288960 5.828930 0.018439\n"
                       "H 2.024414 4.727295 -1.251800\n");

  Utils::AtomCollection structure = Utils::XyzStreamHandler::read(ss);
  calculator.settings().modifyString(Utils::SettingsNames::parameterFilePath, gaff_parameter_file);
  calculator.settings().modifyBool(SwooseUtilities::SettingsNames::detectBondsWithCovalentRadii, true);
  calculator.settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 20.0); // basically no
                                                                                                     // cutoff radius

  /*
   * Add atomic charges file.
   * These charges are taken from OpenBabel in order to reproduce their energy.
   */
  Eigen::VectorXd atomicCharges(8);
  atomicCharges << -0.068144, 0.022715, 0.022715, 0.022715, -0.068144, 0.022715, 0.022715, 0.022715;
  std::ofstream file(atomicChargesFile);
  Utils::matrixToCsv(file, atomicCharges, ',');
  file.close();

  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);
  ASSERT_THROW(calculator.setStructure(structure), std::runtime_error);
  // now set the atomic charges file
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::gaffAtomicChargesFile, atomicChargesFile);
  calculator.setStructure(structure);
  // Clone calculator to test clone interface as well
  auto clonedCalculator = calculator.clone();

  // Calculate and gather results
  Utils::Results results = clonedCalculator->calculate("test calculation ethane");
  double energy = results.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;
  Utils::GradientCollection gradients = results.get<Utils::Property::Gradients>() * Utils::Constants::kCalPerMol_per_hartree;
  Utils::HessianMatrix hessian = results.get<Utils::Property::Hessian>();

  /*
   * Reference taken from OpenBabel implementation.
   * By looking at the verbose output of OpenBabel for each of the interactions,
   * it was concluded that the uncertainty of 0.01 kcal/mol is probably
   * due to rounding errors in OpenBabel.
   */
  ASSERT_THAT(energy, DoubleNear(8.998, 0.01));

  // Numerical gradients
  Utils::PositionCollection positions = structure.getPositions();
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
  for (int k = 0; k < gradients.rows(); ++k) {
    for (int l = 0; l < 3; ++l) {
      ASSERT_THAT(gradients(k, l), DoubleNear(numGrad(k, l), 1e-3));
    }
  }

  // Calculated Hessian fully numerically
  clonedCalculator->setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  Utils::NumericalHessianCalculator numericalHessianCalculator(*clonedCalculator);
  Utils::Results numericalResult = numericalHessianCalculator.calculate();
  const Utils::HessianMatrix& numericalHessian = numericalResult.get<Utils::Property::Hessian>();

  // Assertions
  ASSERT_THAT(hessian.rows(), Eq(numericalHessian.rows()));
  ASSERT_THAT(hessian.cols(), Eq(numericalHessian.cols()));
  for (int k = 0; k < hessian.rows(); ++k) {
    for (int l = 0; l < hessian.cols(); ++l) {
      ASSERT_THAT(hessian(k, l), DoubleNear(numericalHessian(k, l), 1e-5));
    }
  }
}

TEST_F(AGaffMolecularMechanicsTest, GaffCalculationOfAlanineWithReadParameters) {
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;
  calculator.settings().modifyString(Utils::SettingsNames::parameterFilePath, gaff_parameter_file);
  calculator.settings().modifyBool(SwooseUtilities::SettingsNames::detectBondsWithCovalentRadii, true);
  calculator.settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 50.0); // basically no
                                                                                                     // cutoff radius
  calculator.settings().modifyBool(SwooseUtilities::SettingsNames::printContributionsMolecularMechanics, true);
  /*
   * Add atomic charges file with charges taken from OpenBabel.
   */
  Eigen::VectorXd atomicCharges(13);
  atomicCharges << 0.118947, -0.318620, 0.118947, 0.100400, 0.057106, -0.039491, 0.024991, 0.024991, 0.024991, 0.321440,
      -0.249265, -0.479541, 0.295103;
  std::ofstream file(atomicChargesFile);
  Utils::matrixToCsv(file, atomicCharges, ',');
  file.close();
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::gaffAtomicChargesFile, atomicChargesFile);

  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);

  calculator.setStructure(structure);

  // Calculate and gather results
  Utils::Results results = calculator.calculate("test calculation alanine");
  double energy = results.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;
  Utils::GradientCollection gradients = results.get<Utils::Property::Gradients>() * Utils::Constants::kCalPerMol_per_hartree;
  Utils::HessianMatrix hessian = results.get<Utils::Property::Hessian>();

  /*
   * This reference value does not match with the energy obtained from the OpenBabel implementation of GAFF.
   * In OpenBabel, the improper dihedral describing the carboxylic acid group is missing (c3 - o - c - oh).
   * Furthermore, two dihedrals which should be included two potentials (Fourier series) are only included with one
   * by OpenBabel (ho - oh - c - o and hc - c3 - c - o). From the description of the GAFF method, we do not
   * see a reason for omitting these terms.
   * Attempts to get a reasonable reference with the Amber program were not successful yet.
   */
  ASSERT_THAT(energy, DoubleNear(1.97913, 1e-3));

  // Numerical gradients
  Utils::PositionCollection positions = structure.getPositions();
  Utils::GradientCollection numGrad(positions.rows(), 3);

  double stepsize = 1e-5;
  for (int i = 0; i < positions.rows(); ++i) {
    Eigen::RowVector3d pos = positions.row(i);
    Eigen::RowVector3d tmpPos = pos;
    for (int j = 0; j < 3; ++j) {
      pos[j] += stepsize;
      positions.row(i) = pos;

      calculator.setRequiredProperties(Utils::Property::Energy); // This is needed, but why? Bug?
      calculator.modifyPositions(positions);
      results = calculator.calculate("");
      double energyPlus = results.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;

      pos = tmpPos;
      pos[j] -= stepsize;
      positions.row(i) = pos;

      calculator.modifyPositions(positions);
      results = calculator.calculate("");
      double energyMinus = results.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;

      pos = tmpPos;
      positions.row(i) = pos;

      numGrad(i, j) = (energyPlus - energyMinus) / (2.0 * stepsize);
    }
  }

  // Assertions
  for (int k = 0; k < gradients.rows(); ++k) {
    for (int l = 0; l < 3; ++l) {
      ASSERT_THAT(gradients(k, l), DoubleNear(numGrad(k, l), 1e-3));
    }
  }

  // Calculated Hessian fully numerically
  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  Utils::NumericalHessianCalculator numericalHessianCalculator(calculator);
  Utils::Results numericalResult = numericalHessianCalculator.calculate();
  const Utils::HessianMatrix& numericalHessian = numericalResult.get<Utils::Property::Hessian>();

  // Assertions
  ASSERT_THAT(hessian.rows(), Eq(numericalHessian.rows()));
  ASSERT_THAT(hessian.cols(), Eq(numericalHessian.cols()));
  for (int k = 0; k < hessian.rows(); ++k) {
    for (int l = 0; l < hessian.cols(); ++l) {
      ASSERT_THAT(hessian(k, l), DoubleNear(numericalHessian(k, l), 1e-5));
    }
  }
}

TEST_F(AGaffMolecularMechanicsTest, GaffCalculationOfAlanineWithDefaultParameters) {
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;
  calculator.settings().modifyBool(SwooseUtilities::SettingsNames::detectBondsWithCovalentRadii, true);
  calculator.settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 50.0); // basically no
  // cutoff radius
  calculator.settings().modifyBool(SwooseUtilities::SettingsNames::printContributionsMolecularMechanics, true);
  /*
   * Add atomic charges file with charges taken from OpenBabel.
   */
  Eigen::VectorXd atomicCharges(13);
  atomicCharges << 0.118947, -0.318620, 0.118947, 0.100400, 0.057106, -0.039491, 0.024991, 0.024991, 0.024991, 0.321440,
      -0.249265, -0.479541, 0.295103;
  std::ofstream file(atomicChargesFile);
  Utils::matrixToCsv(file, atomicCharges, ',');
  file.close();
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::gaffAtomicChargesFile, atomicChargesFile);

  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian);

  calculator.setStructure(structure);

  // Calculate and gather results
  Utils::Results results = calculator.calculate("test calculation alanine");
  double energy = results.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;
  Utils::GradientCollection gradients = results.get<Utils::Property::Gradients>() * Utils::Constants::kCalPerMol_per_hartree;
  Utils::HessianMatrix hessian = results.get<Utils::Property::Hessian>();

  /*
   * This reference value does not match with the energy obtained from the OpenBabel implementation of GAFF.
   * In OpenBabel, the improper dihedral describing the carboxylic acid group is missing (c3 - o - c - oh).
   * Furthermore, two dihedrals which should be included two potentials (Fourier series) are only included with one
   * by OpenBabel (ho - oh - c - o and hc - c3 - c - o). From the description of the GAFF method, we do not
   * see a reason for omitting these terms.
   * Attempts to get a reasonable reference with the Amber program were not successful yet.
   */
  ASSERT_THAT(energy, DoubleNear(1.97913, 1e-3));

  // Numerical gradients
  Utils::PositionCollection positions = structure.getPositions();
  Utils::GradientCollection numGrad(positions.rows(), 3);

  double stepsize = 1e-5;
  for (int i = 0; i < positions.rows(); ++i) {
    Eigen::RowVector3d pos = positions.row(i);
    Eigen::RowVector3d tmpPos = pos;
    for (int j = 0; j < 3; ++j) {
      pos[j] += stepsize;
      positions.row(i) = pos;

      calculator.setRequiredProperties(Utils::Property::Energy); // This is needed, but why? Bug?
      calculator.modifyPositions(positions);
      results = calculator.calculate("");
      double energyPlus = results.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;

      pos = tmpPos;
      pos[j] -= stepsize;
      positions.row(i) = pos;

      calculator.modifyPositions(positions);
      results = calculator.calculate("");
      double energyMinus = results.get<Utils::Property::Energy>() * Utils::Constants::kCalPerMol_per_hartree;

      pos = tmpPos;
      positions.row(i) = pos;

      numGrad(i, j) = (energyPlus - energyMinus) / (2.0 * stepsize);
    }
  }

  // Assertions
  for (int k = 0; k < gradients.rows(); ++k) {
    for (int l = 0; l < 3; ++l) {
      ASSERT_THAT(gradients(k, l), DoubleNear(numGrad(k, l), 1e-3));
    }
  }

  // Calculated Hessian fully numerically
  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  Utils::NumericalHessianCalculator numericalHessianCalculator(calculator);
  Utils::Results numericalResult = numericalHessianCalculator.calculate();
  const Utils::HessianMatrix& numericalHessian = numericalResult.get<Utils::Property::Hessian>();

  // Assertions
  ASSERT_THAT(hessian.rows(), Eq(numericalHessian.rows()));
  ASSERT_THAT(hessian.cols(), Eq(numericalHessian.cols()));
  for (int k = 0; k < hessian.rows(); ++k) {
    for (int l = 0; l < hessian.cols(); ++l) {
      ASSERT_THAT(hessian(k, l), DoubleNear(numericalHessian(k, l), 1e-5));
    }
  }
}

} // namespace Tests
} // namespace Scine
