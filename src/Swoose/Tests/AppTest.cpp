/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "../App/AppUtils/Tasks.h"
#include "Files/tests_file_location.h"
#include "TestUtilities/MockQmCalculator.h"
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Swoose/MMParametrization/MMParametrizationSettings.h>
#include <Swoose/MMParametrization/Parametrizer.h>
#include <Swoose/MolecularMechanics/SFAM/SfamCalculatorSettings.h>
#include <Swoose/MolecularMechanics/SFAM/SfamMolecularMechanicsCalculator.h>
#include <Swoose/QMMM/QmmmCalculator.h>
#include <Swoose/QMMM/QmmmCalculatorSettings.h>
#include <Utils/GeometryOptimization/QmmmGeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/MolecularDynamics/MolecularDynamics.h>
#include <Utils/MolecularDynamics/MolecularDynamicsSettings.h>
#include <gmock/gmock.h>
#include <yaml-cpp/yaml.h>

using namespace testing;
namespace Scine {
namespace Tests {

/**
 * @class AppTest AppTest.cpp
 * @brief Tests the functionality of task handlers applied in the app. It basically just tests whether
 *        the tasks run through normally without throwing exceptions. The correctness of results is tested
 *        individually for the classes that do the work in other tests.
 * @test
 */
class AppTest : public Test {
 public:
  Core::Log log = Core::Log::silent();
};

TEST_F(AppTest, TaskForSFAMCalculationWorks) {
  MolecularMechanics::SfamMolecularMechanicsCalculator mmCalculator;
  mmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  mmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  mmCalculator.setLog(log);

  Swoose::Tasks::runMMCalculationTask(mmCalculator, melatonin_xyz_file,
                                      Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian, log);
}

TEST_F(AppTest, TaskForSFAMParametrizationWorksWithOrca) {
  MMParametrization::Parametrizer parametrizer;
  auto connFile = Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "Connectivity.dat");
  auto parFile = Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "Parameters.dat");
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, alanine_ref_calc_dir);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, parFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useCsvInputFormat, false);
  parametrizer.setLog(log);

  /*
   * Next line should not cause an error, because the options for reference program
   * should be case-insensitive.
   */
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceProgram, "OrCa");

  Swoose::Tasks::runSFAMParametrizationTask(parametrizer, alanine_xyz_file, log);

  ASSERT_TRUE(boost::filesystem::exists(connFile));
  boost::filesystem::remove(connFile);
  ASSERT_FALSE(boost::filesystem::exists(connFile));
  ASSERT_TRUE(boost::filesystem::exists(parFile));
  boost::filesystem::remove(parFile);
  ASSERT_FALSE(boost::filesystem::exists(parFile));
}

TEST_F(AppTest, TaskForSFAMParametrizationWorksWithTurbomole) {
  MMParametrization::Parametrizer parametrizer;
  auto connFile = Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "Connectivity.dat");
  auto parFile = Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "Parameters.dat");
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, alanine_ref_calc_dir);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, parFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useCsvInputFormat, false);
  parametrizer.setLog(log);

  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceProgram, "Turbomole");

  Swoose::Tasks::runSFAMParametrizationTask(parametrizer, alanine_xyz_file, log);

  ASSERT_TRUE(boost::filesystem::exists(connFile));
  boost::filesystem::remove(connFile);
  ASSERT_FALSE(boost::filesystem::exists(connFile));
  ASSERT_TRUE(boost::filesystem::exists(parFile));
  boost::filesystem::remove(parFile);
  ASSERT_FALSE(boost::filesystem::exists(parFile));
}

TEST_F(AppTest, TaskForClassicalMolecularDynamicsWorks) {
  MolecularMechanics::SfamMolecularMechanicsCalculator mmCalculator;
  mmCalculator.setLog(log);
  mmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  mmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);

  Utils::MolecularDynamics molecularDynamics(mmCalculator);
  molecularDynamics.settings().modifyInt(Utils::SettingsNames::numberOfMDSteps, 9);

  Swoose::Tasks::runMDSimulationTask(molecularDynamics, melatonin_xyz_file, log);

  std::string pathToTrj = "MD_trajectory.xyz";
  std::string pathToEnergies = "MD_energies.dat";
  ASSERT_TRUE(boost::filesystem::exists(pathToTrj));
  boost::filesystem::remove(pathToTrj);
  ASSERT_FALSE(boost::filesystem::exists(pathToTrj));
  ASSERT_TRUE(boost::filesystem::exists(pathToEnergies));

  std::ifstream energiesFile;
  std::string line;
  int nLines;
  energiesFile.open(pathToEnergies);
  for (nLines = 0; std::getline(energiesFile, line); nLines++) {
    continue;
  }
  energiesFile.close();
  ASSERT_THAT(nLines, Eq(11)); // Energies after 9 MD steps + initial energy + comment line at top

  boost::filesystem::remove(pathToEnergies);
  ASSERT_FALSE(boost::filesystem::exists(pathToEnergies));
}

TEST_F(AppTest, TaskForHybridMolecularDynamicsWorks) {
  auto& manager = Core::ModuleManager::getInstance();

  if (!manager.moduleLoaded("MockModule")) {
    manager.load("Test/TestUtilities");
  }

  auto mmCalculator = std::make_shared<MolecularMechanics::SfamMolecularMechanicsCalculator>();
  auto mmCalculatorAsCalculator = std::dynamic_pointer_cast<Core::Calculator>(mmCalculator);
  auto qmCalculator = manager.get<Core::Calculator>("MOCK-QM", "MockModule");
  Qmmm::QmmmCalculator qmmmCalculator;
  qmmmCalculator.setUnderlyingCalculators(qmCalculator, mmCalculatorAsCalculator);
  qmmmCalculator.setLog(log);
  qmmmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  qmmmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  qmmmCalculator.settings().modifyBool(SwooseUtilities::SettingsNames::electrostaticEmbedding, false);
  const std::vector<int> qmAtomsList{{30, 31, 32, 23, 22, 29, 21, 28, 20, 26, 27, 19, 24, 25}};
  qmmmCalculator.settings().modifyIntList(SwooseUtilities::SettingsNames::qmAtomsList, qmAtomsList);

  Utils::MolecularDynamics molecularDynamics(qmmmCalculator);
  molecularDynamics.settings().modifyInt(Utils::SettingsNames::numberOfMDSteps, 7);

  Swoose::Tasks::runMDSimulationTask(molecularDynamics, melatonin_xyz_file, log);

  std::string pathToTrj = "MD_trajectory.xyz";
  std::string pathToEnergies = "MD_energies.dat";
  ASSERT_TRUE(boost::filesystem::exists(pathToTrj));
  boost::filesystem::remove(pathToTrj);
  ASSERT_FALSE(boost::filesystem::exists(pathToTrj));
  ASSERT_TRUE(boost::filesystem::exists(pathToEnergies));

  std::ifstream energiesFile;
  std::string line;
  int nLines;
  energiesFile.open(pathToEnergies);
  for (nLines = 0; std::getline(energiesFile, line); nLines++) {
    continue;
  }
  energiesFile.close();
  ASSERT_THAT(nLines, Eq(9)); // Energies after 7 MD steps + initial energy + comment line at top

  boost::filesystem::remove(pathToEnergies);
  ASSERT_FALSE(boost::filesystem::exists(pathToEnergies));
}

TEST_F(AppTest, TaskForQmmmWorks) {
  auto& manager = Core::ModuleManager::getInstance();

  if (!manager.moduleLoaded("MockModule")) {
    manager.load("Test/TestUtilities");
  }

  auto mmCalculator = std::make_shared<MolecularMechanics::SfamMolecularMechanicsCalculator>();
  auto mmCalculatorAsCalculator = std::dynamic_pointer_cast<Core::Calculator>(mmCalculator);
  auto qmCalculator = manager.get<Core::Calculator>("MOCK-QM", "MockModule");
  Qmmm::QmmmCalculator qmmmCalculator;
  qmmmCalculator.setUnderlyingCalculators(qmCalculator, mmCalculatorAsCalculator);
  qmmmCalculator.setLog(log);
  qmmmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  qmmmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  qmmmCalculator.settings().modifyBool(SwooseUtilities::SettingsNames::electrostaticEmbedding, false);

  Utils::AtomCollection testStructure = Utils::ChemicalFileHandler::read(melatonin_xyz_file).first;
  qmmmCalculator.setStructure(testStructure);

  const std::vector<int> qmAtoms{{30, 31, 32, 23, 22, 29, 21, 28, 20, 26, 27, 19, 24, 25}};
  YAML::Node yamlNode;
  yamlNode["qm_atoms"] = qmAtoms;
  Swoose::Tasks::runQmmmCalculationTask(qmmmCalculator, melatonin_xyz_file,
                                        Utils::Property::Energy | Utils::Property::Gradients, log, yamlNode);
}

TEST_F(AppTest, TaskForSFAMOptimizationWorks) {
  MolecularMechanics::SfamMolecularMechanicsCalculator mmCalculator;
  mmCalculator.setLog(log);
  mmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, water_param_file);
  mmCalculator.settings().modifyBool(SwooseUtilities::SettingsNames::detectBondsWithCovalentRadii, true);

  YAML::Node yamlNode;
  auto optimizer = std::make_unique<Utils::GeometryOptimizer<Utils::Bfgs>>(mmCalculator);
  optimizer->check.maxIter = 100;
  optimizer->optimizer.useTrustRadius = true;
  Swoose::Tasks::runMMOptimizationTask(mmCalculator, *optimizer, water_xyz_file, log, yamlNode);

  std::string pathToOptResults = "opt_results";
  ASSERT_TRUE(boost::filesystem::exists(pathToOptResults));
  boost::filesystem::remove_all(pathToOptResults);
  ASSERT_FALSE(boost::filesystem::exists(pathToOptResults));
}

TEST_F(AppTest, TaskForQmmmOptimizationWorks) {
  auto& manager = Core::ModuleManager::getInstance();

  if (!manager.moduleLoaded("MockModule")) {
    manager.load("Test/TestUtilities");
  }

  auto mmCalculator = std::make_shared<MolecularMechanics::SfamMolecularMechanicsCalculator>();
  auto mmCalculatorAsCalculator = std::dynamic_pointer_cast<Core::Calculator>(mmCalculator);
  auto qmCalculator = manager.get<Core::Calculator>("MOCK-QM", "MockModule");
  Qmmm::QmmmCalculator qmmmCalculator;
  qmmmCalculator.setUnderlyingCalculators(qmCalculator, mmCalculatorAsCalculator);
  qmmmCalculator.setLog(log);
  qmmmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, melatonin_param_file);
  qmmmCalculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, melatonin_connectivity_file);
  qmmmCalculator.settings().modifyBool(SwooseUtilities::SettingsNames::electrostaticEmbedding, false);

  const std::vector<int> qmAtoms{{30, 31, 32, 23, 22, 29, 21, 28, 20, 26, 27, 19, 24, 25}};
  YAML::Node yamlNode;
  yamlNode["qm_atoms"] = qmAtoms;

  auto optimizer = std::make_unique<Utils::QmmmGeometryOptimizer<Utils::Bfgs>>(qmmmCalculator);

  optimizer->maxMacrocycles = 2;
  optimizer->maxFullOptMicrocycles = 3;
  optimizer->maxEnvOptMicrocycles = 5;
  optimizer->coordinateSystem = Utils::CoordinateSystem::Cartesian;

  // With the mock QM calculator, it won't converge.
  std::string exceptionString = "";
  try {
    Swoose::Tasks::runQmmmOptimizationTask<Utils::Bfgs>(qmmmCalculator, *optimizer, melatonin_xyz_file, log, yamlNode);
  }
  catch (const std::runtime_error& e) {
    exceptionString = e.what();
  }

  ASSERT_STREQ(exceptionString.c_str(),
               "The structure optimization did not converge within the maximum number of cycles.");

  // Directory should still exist.
  std::string pathToOptResults = "opt_results";
  ASSERT_TRUE(boost::filesystem::exists(pathToOptResults));
  boost::filesystem::remove_all(pathToOptResults);
  ASSERT_FALSE(boost::filesystem::exists(pathToOptResults));
}

} // namespace Tests
} // namespace Scine
