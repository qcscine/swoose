/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include <Swoose/MolecularMechanics/GAFF/GaffCalculatorSettings.h>
#include <Swoose/MolecularMechanics/GAFF/GaffMolecularMechanicsCalculator.h>
#include <Swoose/MolecularMechanics/SFAM/SfamCalculatorSettings.h>
#include <Swoose/MolecularMechanics/SFAM/SfamMolecularMechanicsCalculator.h>
#include <Swoose/QMMM/InteractionTermEliminator.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>
#include <numeric>

using namespace testing;
namespace Scine {
using namespace Qmmm;
using namespace MolecularMechanics;
namespace Tests {

/**
 * @class SfamInteractionTermEliminatorTest InteractionTermEliminatorTest.cpp
 * @brief Tests the functionality of the InteractionTermEliminator class applied in QM/MM with SFAM as MM.
 * @test
 */
class SfamInteractionTermEliminatorTest : public Test {
 public:
  std::shared_ptr<SfamMolecularMechanicsCalculator> calculator;
  Core::Log silentLogger = Core::Log::silent();
  void SetUp() override {
    calculator = std::make_shared<SfamMolecularMechanicsCalculator>();
    calculator->setLog(silentLogger);
  }
};

/**
 * @class GaffInteractionTermEliminatorTest InteractionTermEliminatorTest.cpp
 * @brief Tests the functionality of the InteractionTermEliminator class applied in QM/MM with GAFF as MM.
 * @test
 */
class GaffInteractionTermEliminatorTest : public Test {
 public:
  std::shared_ptr<GaffMolecularMechanicsCalculator> calculator;
  Core::Log silentLogger = Core::Log::silent();
  std::string atomicChargesFile = "atomic_charges_for_gaff_test.dat";
  void SetUp() override {
    calculator = std::make_shared<GaffMolecularMechanicsCalculator>();
    calculator->setLog(silentLogger);
  }
  void TearDown() override {
    boost::filesystem::remove_all(atomicChargesFile);
  }
};

TEST_F(SfamInteractionTermEliminatorTest, SfamInteractionTermsAreCompletelyDisabledForMelatonin) {
  Utils::AtomCollection testStructure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;

  calculator->settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  calculator->settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  calculator->settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 20.0); // basically no
                                                                                                      // cutoff radius

  calculator->setStructure(testStructure);

  auto results = calculator->calculate("test calculation 1");
  double initialEnergy = results.get<Utils::Property::Energy>();

  // When no interaction terms are eliminated, then the absolute value of the energy should be larger than zero.
  ASSERT_TRUE(std::abs(initialEnergy) > 1e-3);

  std::vector<int> qmAtoms(testStructure.size());
  std::iota(qmAtoms.begin(), qmAtoms.end(), 0);

  // Next line not needed (see test below), but it should work, too.
  auto baseCalcPtr = std::static_pointer_cast<MolecularMechanicsCalculator>(calculator);

  InteractionTermEliminator completeEliminator(qmAtoms, baseCalcPtr);
  completeEliminator.eliminateInteractionTerms(true);

  // Calculation with completely disabling all interaction terms.
  results = calculator->calculate("test calculation 2");
  ASSERT_THAT(results.get<Utils::Property::Energy>(), DoubleNear(0.0, 1e-6));

  completeEliminator.reset();
  results = calculator->calculate("test calculation 3");

  // After reverting to the initial state, the energy should equal the initial energy.
  ASSERT_THAT(results.get<Utils::Property::Energy>(), DoubleNear(initialEnergy, 1e-6));

  completeEliminator.eliminateInteractionTerms(false);
  results = calculator->calculate("test calculation 4");
  ASSERT_THAT(results.get<Utils::Property::Energy>(), DoubleNear(0.0, 1e-6));
}

TEST_F(SfamInteractionTermEliminatorTest, SfamInteractionTermsArePartiallyDisabledForMelatonin) {
  Utils::AtomCollection testStructure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;

  calculator->settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  calculator->settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  calculator->settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 20.0); // basically no
                                                                                                      // cutoff radius

  calculator->setStructure(testStructure);

  std::vector<int> qmAtoms = {30, 31, 32, 23, 22, 29, 21, 28, 20, 26, 27, 19, 24, 25};

  InteractionTermEliminator partialEliminator(qmAtoms, calculator);
  partialEliminator.eliminateInteractionTerms(true); // electrostatic embedding: on

  auto results = calculator->calculate("test calculation 1");
  auto energy1 = results.get<Utils::Property::Energy>();
  ASSERT_TRUE(std::abs(energy1 * Utils::Constants::kCalPerMol_per_hartree) > 1.0);

  partialEliminator.reset();
  partialEliminator.eliminateInteractionTerms(false); // electrostatic embedding: off

  results = calculator->calculate("test calculation 2");
  auto energy2 = results.get<Utils::Property::Energy>();
  ASSERT_TRUE(std::abs(energy2 * Utils::Constants::kCalPerMol_per_hartree) > 1.0);

  ASSERT_TRUE(std::abs(energy2 - energy1) * Utils::Constants::kCalPerMol_per_hartree > 0.5);
}

TEST_F(GaffInteractionTermEliminatorTest, GaffInteractionTermsAreCompletelyDisabledForAlanine) {
  Utils::AtomCollection testStructure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;

  calculator->settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, gaff_parameter_file);
  calculator->settings().modifyBool(SwooseUtilities::SettingsNames::detectBondsWithCovalentRadii, true);
  calculator->settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 20.0); // basically no
                                                                                                      // cutoff radius

  calculator->setStructure(testStructure);

  auto results = calculator->calculate("test calculation 1");
  double initialEnergy = results.get<Utils::Property::Energy>();

  // When no interaction terms are eliminated, then the absolute value of the energy should be larger than zero.
  ASSERT_TRUE(std::abs(initialEnergy) > 1e-3);

  std::vector<int> qmAtoms(testStructure.size());
  std::iota(qmAtoms.begin(), qmAtoms.end(), 0);

  // Next line not needed (see test below), but it should work, too.
  auto baseCalcPtr = std::static_pointer_cast<MolecularMechanicsCalculator>(calculator);

  InteractionTermEliminator completeEliminator(qmAtoms, baseCalcPtr);
  completeEliminator.eliminateInteractionTerms(true);

  // Calculation with completely disabling all interaction terms.
  results = calculator->calculate("test calculation 2");
  ASSERT_THAT(results.get<Utils::Property::Energy>(), DoubleNear(0.0, 1e-6));

  completeEliminator.reset();
  results = calculator->calculate("test calculation 3");

  // After reverting to the initial state, the energy should equal the initial energy.
  ASSERT_THAT(results.get<Utils::Property::Energy>(), DoubleNear(initialEnergy, 1e-6));

  completeEliminator.eliminateInteractionTerms(false);
  results = calculator->calculate("test calculation 4");
  ASSERT_THAT(results.get<Utils::Property::Energy>(), DoubleNear(0.0, 1e-6));
}

TEST_F(GaffInteractionTermEliminatorTest, GaffInteractionTermsArePartiallyDisabledForAlanine) {
  Utils::AtomCollection testStructure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;

  calculator->settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, gaff_parameter_file);
  calculator->settings().modifyBool(SwooseUtilities::SettingsNames::detectBondsWithCovalentRadii, true);
  calculator->settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 20.0); // basically no
                                                                                                      // cutoff radius

  /*
   * Add atomic charges file.
   * These charges are taken from OpenBabel.
   */
  Eigen::VectorXd atomicCharges(13);
  atomicCharges << 0.118947, -0.318620, 0.118947, 0.100400, 0.057106, -0.039491, 0.024991, 0.024991, 0.024991, 0.321440,
      -0.249265, -0.479541, 0.295103;
  /*
   * Artificially double charges such that the effect measured below is more pronounced
   * and can be measured against the same thresholds as in the SFAM test example.
   */
  atomicCharges *= 2.0;
  std::ofstream file(atomicChargesFile);
  Utils::matrixToCsv(file, atomicCharges, ',');
  file.close();
  calculator->settings().modifyString(SwooseUtilities::SettingsNames::gaffAtomicChargesFile, atomicChargesFile);

  calculator->setStructure(testStructure);

  std::vector<int> qmAtoms = {3, 4, 5, 6, 7, 8};

  InteractionTermEliminator partialEliminator(qmAtoms, calculator);
  partialEliminator.eliminateInteractionTerms(true); // electrostatic embedding: on

  auto results = calculator->calculate("test calculation 1");
  auto energy1 = results.get<Utils::Property::Energy>();
  ASSERT_TRUE(std::abs(energy1 * Utils::Constants::kCalPerMol_per_hartree) > 1.0);

  partialEliminator.reset();
  partialEliminator.eliminateInteractionTerms(false); // electrostatic embedding: off

  results = calculator->calculate("test calculation 2");
  auto energy2 = results.get<Utils::Property::Energy>();
  ASSERT_TRUE(std::abs(energy2 * Utils::Constants::kCalPerMol_per_hartree) > 1.0);

  ASSERT_TRUE(std::abs(energy2 - energy1) * Utils::Constants::kCalPerMol_per_hartree > 0.5);
}

} // namespace Tests
} // namespace Scine
