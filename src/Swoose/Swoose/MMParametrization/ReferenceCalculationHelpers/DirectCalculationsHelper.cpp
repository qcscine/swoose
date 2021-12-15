/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DirectCalculationsHelper.h"
#include "../MMParametrizationSettings.h"
#include "../ParametrizationData.h"
#include "BasicJobSubmissionHelper.h"
#include <Core/ModuleManager.h>
#include <Utils/ExternalQC/Gaussian/GaussianCalculator.h>
#include <Utils/ExternalQC/Gaussian/GaussianCalculatorSettings.h>
#include <Utils/ExternalQC/Orca/OrcaCalculator.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/Yaml.h>
#include <Utils/Properties/AtomicCharges/ChargeModel5.h>
#include <yaml-cpp/yaml.h>
#include <boost/filesystem.hpp>

namespace Scine {
namespace MMParametrization {
namespace DirectCalculationsHelper {

namespace {
static constexpr const char* optimizationTrajectoryFile = "optimization_trajectory.xyz";
} // namespace

void performReferenceCalculations(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings,
                                  std::string baseWorkingDir, Core::Log& log) {
  if (data.vectorOfStructures.size() != 1)
    throw std::runtime_error("Cannot perform direct calculations for fragments.");

  // Get a module manager
  auto& manager = Core::ModuleManager::getInstance();

  // Get the calculator settings:
  auto yamlSettingsPath = settings->getString(SwooseUtilities::SettingsNames::yamlSettingsFilePath);
  YAML::Node yamlSettings;
  if (!yamlSettingsPath.empty())
    yamlSettings = YAML::LoadFile(yamlSettingsPath);

  // Set up the calculator:
  auto referenceMethod = settings->getString(SwooseUtilities::SettingsNames::referenceMethod);
  auto referenceBasisSet = settings->getString(SwooseUtilities::SettingsNames::referenceBasisSet);

  auto referenceProgram = settings->getString(SwooseUtilities::SettingsNames::referenceProgram);
  // For the next line, referenceProgram will be all-lowercase
  auto methodFamily = BasicJobSubmissionHelper::determineMethodFamily(referenceMethod, referenceProgram);

  // Convert capitalization
  referenceProgram.at(0) = std::toupper(referenceProgram.at(0));
  std::transform(methodFamily.begin(), methodFamily.end(), methodFamily.begin(), ::toupper);

  std::shared_ptr<Core::Calculator> calculator;
  try {
    calculator = manager.get<Core::Calculator>(Core::Calculator::supports(methodFamily), referenceProgram);
  }
  catch (const std::runtime_error& error) {
    throw std::runtime_error(
        "The reference calculator could not be loaded via the module system.\nCheck: (i) whether the requested "
        "calculator is available (see manual), (ii) that you have installed all "
        "relevant modules, and (iii) that all necessary environment variables are set (e.g., ORCA_BINARY_PATH or "
        "TURBODIR).");
  }

  // Set logger
  Core::Log warningLog = Core::Log::silent();
  warningLog.warning.add("cerr", Core::Log::cerrSink());
  warningLog.error.add("cerr", Core::Log::cerrSink());
  calculator->setLog(warningLog);

  Utils::nodeToSettings(calculator->settings(), yamlSettings, true);
  calculator->settings().modifyInt(Utils::SettingsNames::molecularCharge, data.vectorOfChargesAndMultiplicities.at(0).first);
  calculator->settings().modifyInt(Utils::SettingsNames::spinMultiplicity, data.vectorOfChargesAndMultiplicities.at(0).second);
  if (calculator->settings().valueExists(Utils::ExternalQC::SettingsNames::baseWorkingDirectory))
    calculator->settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, baseWorkingDir);
  if (calculator->settings().valueExists(Utils::SettingsNames::method))
    calculator->settings().modifyString(Utils::SettingsNames::method, referenceMethod);
  if (calculator->settings().valueExists(Utils::SettingsNames::basisSet))
    calculator->settings().modifyString(Utils::SettingsNames::basisSet, referenceBasisSet);

  Utils::GeometryOptimizer<Utils::Bfgs> optimizer(*calculator);
  Utils::AtomCollection structure = *data.vectorOfStructures.at(0);

  // Trajectory for optimized structure
  std::ofstream trajectory(optimizationTrajectoryFile, std::ofstream::out);
  Utils::XyzStreamHandler writer;

  // Variable for the observer
  double oldEnergy = 0.0;
  // Observer function
  auto obsFunc = [&](const int& cycle, const double& energy, const Eigen::VectorXd& /* params */) {
    if (cycle == 1) {
      log.debug.printf("%7s %16s %16s\n", "Cycle", "Energy", "Energy Diff.");
    }
    log.debug.printf("%7d %+16.9f %+16.9f\n", cycle, energy, energy - oldEnergy);
    oldEnergy = energy;
    auto currStructure = calculator->getStructure();
    writer.write(trajectory, *currStructure);
  };

  // Add observer
  optimizer.addObserver(obsFunc);

  auto optimizerSettings = optimizer.getSettings();
  Utils::nodeToSettings(optimizerSettings, yamlSettings, true);
  optimizer.setSettings(optimizerSettings);

  log.output << "Optimizing structure...";
  log.debug << Core::Log::nl;
  log.output << Core::Log::endl;

  int cycles = 0;
  try {
    cycles = optimizer.optimize(structure, log);
  }
  catch (const std::exception& e) {
    trajectory.close();
    throw std::runtime_error(std::string(e.what()) + "\n" + "The trajectory of the optimization can be found here: " +
                             std::string(optimizationTrajectoryFile));
  }

  if (cycles == optimizer.check.maxIter) {
    trajectory.close();
    throw std::runtime_error("The structure optimization did not converge within the maximum number of cycles.\nThe "
                             "trajectory of the optimization can be found here: " +
                             std::string(optimizationTrajectoryFile));
  }

  // Remove optimization trajectory file if no error occurred:
  boost::filesystem::remove_all(optimizationTrajectoryFile);

  data.vectorOfOptimizedStructures.resize(1);
  data.vectorOfOptimizedStructures[0] = std::make_unique<Utils::AtomCollection>(structure);
  log.output << "Done." << Core::Log::nl;
  log.debug << Core::Log::nl;

  // Perform Hessian calculation
  log.output << "Calculating Hessian..." << Core::Log::endl;
  calculator->setRequiredProperties(Utils::Property::Hessian);
  calculator->modifyPositions(data.vectorOfOptimizedStructures.at(0)->getPositions());
  const auto& hessianResults = calculator->calculate("");
  data.vectorOfHessians[0] = std::make_unique<Utils::HessianMatrix>(hessianResults.get<Utils::Property::Hessian>());
  log.output << "Done." << Core::Log::nl;

  // Perform Mayer bond orders calculation if requested
  if (settings->getBool(SwooseUtilities::SettingsNames::refineConnectivity)) {
    log.output << "Calculating bond orders..." << Core::Log::endl;
    calculator->modifyPositions(data.vectorOfStructures.at(0)->getPositions());
    calculator->setRequiredProperties(Utils::Property::Energy | Utils::Property::BondOrderMatrix);
    const auto& bondOrdersResults = calculator->calculate("");
    data.vectorOfBondOrderCollections.resize(1);
    data.vectorOfBondOrderCollections[0] =
        std::make_unique<Utils::BondOrderCollection>(bondOrdersResults.get<Utils::Property::BondOrderMatrix>());
    log.output << "Done." << Core::Log::endl;
  }

  if (!settings->getBool(SwooseUtilities::SettingsNames::useGaussianOptionKey)) {
    // Perform atomic charges calculation, Hirshfeld if ORCA is the program
    log.output << "Calculating atomic charges (Hirshfeld if ORCA is selected program, Loewdin if TURBOMOLE is selected "
                  "program, otherwise Mulliken)..."
               << Core::Log::endl;
    calculator->modifyPositions(data.vectorOfOptimizedStructures.at(0)->getPositions());
    calculator->setRequiredProperties(Utils::Property::Energy | Utils::Property::AtomicCharges);
    const auto& chargesResults = calculator->calculate("");
    data.atomicChargesForEachFragment.resize(1);
    std::vector<double> atomicCharges = chargesResults.get<Utils::Property::AtomicCharges>();

    // Apply CM5 correction if requested by setting
    if (settings->getBool(SwooseUtilities::SettingsNames::convertChargesCm5)) {
      log.output << "Applying CM5 correction..." << Core::Log::endl;
      data.atomicChargesForEachFragment[0] =
          Utils::ChargeModel5::calculateCm5Charges(atomicCharges, *data.vectorOfOptimizedStructures.at(0));
    }
    else {
      data.atomicChargesForEachFragment[0] = atomicCharges;
    }
    log.output << "Done." << Core::Log::endl;
  }
}

void calculateAtomicChargesWithGaussian(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings,
                                        std::string baseWorkingDir, Core::Log& log) {
  if (!settings->getBool(SwooseUtilities::SettingsNames::useGaussianOptionKey))
    return; // This function should not be called in this case, but double-check to be sure.

  if (data.vectorOfStructures.size() != 1)
    throw std::runtime_error("Cannot calculate atomic charges directly for fragments.");

  // Get the calculator settings:
  auto yamlSettingsPath = settings->getString(SwooseUtilities::SettingsNames::yamlSettingsFilePath);
  YAML::Node yamlSettings;
  if (!yamlSettingsPath.empty())
    yamlSettings = YAML::LoadFile(yamlSettingsPath);

  log.output << "Calculating CM5 charges directly with Gaussian..." << Core::Log::endl;
  Utils::ExternalQC::GaussianCalculator gaussianCalculator;

  // Set the correct settings for the Gaussian calculation
  Utils::nodeToSettings(gaussianCalculator.settings(), yamlSettings, true);
  auto gaussianMethod = settings->getString(SwooseUtilities::SettingsNames::gaussianMethod);
  auto gaussianBasisSet = settings->getString(SwooseUtilities::SettingsNames::gaussianBasisSet);
  gaussianCalculator.settings().modifyString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, baseWorkingDir);
  gaussianCalculator.settings().modifyString(Utils::SettingsNames::method, gaussianMethod);
  gaussianCalculator.settings().modifyString(Utils::SettingsNames::basisSet, gaussianBasisSet);
  gaussianCalculator.settings().modifyInt(Utils::SettingsNames::molecularCharge,
                                          data.vectorOfChargesAndMultiplicities.at(0).first);
  gaussianCalculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity,
                                          data.vectorOfChargesAndMultiplicities.at(0).second);

  // Perform calculation
  gaussianCalculator.setStructure(*data.vectorOfOptimizedStructures.at(0));
  gaussianCalculator.setRequiredProperties(Utils::Property::AtomicCharges);
  const auto& results = gaussianCalculator.calculate("");
  data.atomicChargesForEachFragment.resize(1);
  data.atomicChargesForEachFragment[0] = results.get<Utils::Property::AtomicCharges>();
  log.output << "Done." << Core::Log::endl;
}

} // namespace DirectCalculationsHelper
} // namespace MMParametrization
} // namespace Scine
