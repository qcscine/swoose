/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include "TestUtilities/MockQmCalculator.h"
#include <Core/ModuleManager.h>
#include <Swoose/MolecularMechanics/SFAM/SfamCalculatorSettings.h>
#include <Swoose/MolecularMechanics/SFAM/SfamMolecularMechanicsCalculator.h>
#include <Swoose/QMMM/QmmmCalculator.h>
#include <Swoose/QMMM/QmmmCalculatorSettings.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <gmock/gmock.h>
#include <boost/filesystem/operations.hpp>

using namespace testing;
namespace Scine {
using namespace Qmmm;
namespace Tests {

/**
 * @class QmmmTest QmmmTest.cpp
 * @brief Tests the functionality of some helper functions which are applied in the QM/MM calculator.
 * @test
 */
class QmmmTest : public Test {
 public:
  QmmmCalculator calculator;
  static constexpr const char* qmRegionXyzFile = "created_qm_region_for_test.xyz";
  Utils::AtomCollection testStructure;
  std::vector<int> qmAtoms;

  void SetUp() override {
    testStructure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;
    qmAtoms = {30, 31, 32, 23, 22, 29, 21, 28, 20, 26, 27, 19, 24, 25};

    auto& manager = Core::ModuleManager::getInstance();

    if (!manager.moduleLoaded("MockModule")) {
      manager.load("Test/TestUtilities");
    }

    auto mmCalculator = std::make_shared<MolecularMechanics::SfamMolecularMechanicsCalculator>();
    auto mmCalculatorAsCalculator = std::dynamic_pointer_cast<Core::Calculator>(mmCalculator);
    auto qmCalculator = manager.get<Core::Calculator>("MOCK-QM", "MockModule");
    calculator.setUnderlyingCalculators({qmCalculator, mmCalculatorAsCalculator});
  }

  ~QmmmTest() override {
    boost::filesystem::remove(qmRegionXyzFile);
  }
};

TEST_F(QmmmTest, QmmmSettingsWorkCorrectly) {
  using namespace SwooseUtilities;
  calculator.settings().modifyIntList(Utils::SettingsNames::qmAtomsList, std::vector<int>{{1, 2, 3}});
  calculator.settings().modifyString(SettingsNames::chargeRedistributionKey, OptionNames::redistributedChargeOption);
  calculator.settings().modifyString(Utils::SettingsNames::parameterFilePath, "params.dat");
  calculator.settings().modifyString(SettingsNames::connectivityFilePath, "conn.dat");
  calculator.settings().modifyString(SettingsNames::qmRegionXyzFile, qmRegionXyzFile);
  calculator.settings().modifyBool(SettingsNames::hydrogenBondCorrection, false);
  calculator.settings().modifyBool(Utils::SettingsNames::electrostaticEmbedding, false);
  calculator.settings().modifyDouble(SettingsNames::nonCovalentCutoffRadius, 12.5);
  calculator.settings().modifyBool(SettingsNames::detectBondsWithCovalentRadii, false);
  calculator.settings().modifyBool(SettingsNames::printContributionsMolecularMechanics, true);
  calculator.settings().modifyBool(SettingsNames::onlyCalculateBondedContribution, false);

  ASSERT_THAT(calculator.settings().getIntList(Utils::SettingsNames::qmAtomsList), Eq(std::vector<int>{{1, 2, 3}}));
  ASSERT_THAT(calculator.settings().getString(SettingsNames::chargeRedistributionKey),
              Eq(OptionNames::redistributedChargeOption));
  ASSERT_THAT(calculator.settings().getString(Utils::SettingsNames::parameterFilePath), Eq("params.dat"));
  ASSERT_THAT(calculator.settings().getString(SettingsNames::connectivityFilePath), Eq("conn.dat"));
  ASSERT_THAT(calculator.settings().getString(SettingsNames::qmRegionXyzFile), Eq(qmRegionXyzFile));
  ASSERT_FALSE(calculator.settings().getBool(SettingsNames::hydrogenBondCorrection));
  ASSERT_FALSE(calculator.settings().getBool(Utils::SettingsNames::electrostaticEmbedding));
  ASSERT_THAT(calculator.settings().getDouble(SettingsNames::nonCovalentCutoffRadius), DoubleNear(12.5, 1e-6));
  ASSERT_FALSE(calculator.settings().getBool(SettingsNames::detectBondsWithCovalentRadii));
  ASSERT_TRUE(calculator.settings().getBool(SettingsNames::printContributionsMolecularMechanics));
  ASSERT_FALSE(calculator.settings().getBool(SettingsNames::onlyCalculateBondedContribution));
}

TEST_F(QmmmTest, QmmmCalculationWorksForMelatonin) {
  calculator.settings().modifyString(Utils::SettingsNames::parameterFilePath, melatonin_param_file);
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  calculator.settings().modifyBool(SwooseUtilities::SettingsNames::calculateReducedQmMmEnergy, false);
  calculator.settings().modifyString(SwooseUtilities::SettingsNames::qmRegionXyzFile, qmRegionXyzFile);
  calculator.settings().modifyIntList(Utils::SettingsNames::qmAtomsList, qmAtoms);

  auto log = Core::Log::silent();
  calculator.setLog(log);

  calculator.setStructure(testStructure);
  calculator.setStructure(testStructure); // a second setStructure() should not break anything

  calculator.settings().modifyInt(Swoose::SettingsNames::mockIntSettingKey, 12);
  calculator.settings().modifyString(Swoose::SettingsNames::mockStringSettingKey, "updated string setting");

  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  auto calculators = calculator.getUnderlyingCalculators();
  // calculators are set correctly
  ASSERT_TRUE(std::find(calculators.begin(), calculators.end(), nullptr) == calculators.end());

  auto results = calculator.calculate("test calculation");
  auto energy = results.get<Utils::Property::Energy>();
  Utils::GradientCollection gradients = results.get<Utils::Property::Gradients>();

  ASSERT_THAT(energy, DoubleNear(200.0, 0.1));
  ASSERT_TRUE(std::abs(200.0 - energy) > 0.01);

  for (int i = 0; i < gradients.rows(); ++i) {
    if (i != 2 && i != 19) { // The gradients of these two atoms are coupled
      for (int j = 0; j < 3; ++j) {
        if (std::find(qmAtoms.begin(), qmAtoms.end(), i) != qmAtoms.end()) {
          ASSERT_THAT(gradients(i, j), DoubleNear(100.0, 0.1));
        }
        else {
          ASSERT_THAT(gradients(i, j), DoubleNear(0.0, 0.1));
        }
      }
    }
    else {
      for (int j = 0; j < 3; ++j) {
        ASSERT_FALSE(std::abs(gradients(i, j)) < 0.1);
        ASSERT_FALSE(std::abs(gradients(i, j) - 100.0) < 0.1);
      }
    }
  }

  ASSERT_THAT(results.get<Utils::Property::Description>(), Eq("test calculation"));
  ASSERT_TRUE(boost::filesystem::exists(qmRegionXyzFile));
  ASSERT_FALSE(boost::filesystem::exists("environment_pointcharges.pc"));
  // reduced QM/MM energy is correctly stored in the results
  calculator.settings().modifyBool(SwooseUtilities::SettingsNames::calculateReducedQmMmEnergy, true);
  results = calculator.calculate("test calculation");
  auto reduced_energy = results.get<Utils::Property::Energy>();

  ASSERT_THAT(reduced_energy, DoubleNear(200.0, 0.01));
  ASSERT_TRUE(std::abs(200.0 - reduced_energy) > 0.001);

  // calculate partial hessian
  calculator.setRequiredProperties(Utils::Property::Energy | Utils::Property::PartialHessian);
  auto results2 = calculator.calculate("test calculation");
  Utils::PartialHessian partialHessian = results2.get<Utils::Property::PartialHessian>();
  int hessianSize = qmAtoms.size() * 3;
  ASSERT_THAT(partialHessian.getMatrix().rows(), Eq(hessianSize));
  ASSERT_THAT(partialHessian.getMatrix().cols(), Eq(hessianSize));

  auto partialHessianIndices = partialHessian.getIndices();
  for (long unsigned int i = 0; i < partialHessianIndices.size(); i++) {
    ASSERT_THAT(partialHessianIndices.at(i), Eq(qmAtoms.at(i)));
  }
}

} // namespace Tests
} // namespace Scine
