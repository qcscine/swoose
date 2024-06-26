/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_TASKMANAGEMENT_H
#define SWOOSE_TASKMANAGEMENT_H

#include "Tasks.h"
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Swoose/QMMM/QmmmCalculator.h>
#include <Swoose/StructurePreparation/StructureProcessor.h>
#include <Swoose/Utilities/CalculatorOptions.h>
#include <Swoose/Utilities/SettingsNames.h>
#include <Utils/GeometryOptimization/QmmmGeometryOptimizer.h>
#include <Utils/IO/Yaml.h>
#include <Utils/MolecularDynamics/MolecularDynamics.h>
#include <yaml-cpp/yaml.h>

namespace Scine {
namespace Swoose {
namespace TaskManagement {

/**
 * @brief Helper function to load MM calculator via module system with correct error message.
 * @param mmModel The model name of the calculator.
 * @param manager The module manager.
 * @return The MM calculator.
 */
std::shared_ptr<Core::Calculator> getMmCalculatorFromModel(std::string mmModel, Core::ModuleManager& manager) {
  std::shared_ptr<Core::Calculator> calculator;
  try {
    calculator = manager.get<Core::Calculator>(mmModel);
  }
  catch (const std::runtime_error& e) {
    throw std::runtime_error("MM model could not be loaded via the module system.");
  }
  return calculator;
}

/**
 * @brief Manages all tasks of the app.
 * @param manager The module manager.
 * @param mode The mode of the app: parametrize, calculate, md, optimize.
 * @param quantum Whether to use QM/MM instead of only MM.
 * @param hessianRequired Whether the calculation of the Hessian matrix is required.
 * @param structureFile The path to the XYZ file defining the molecular structure.
 * @param yamlNode The settings as a YAML node.
 * @param log The logger.
 */
void manageTasks(Core::ModuleManager& manager, std::string mode, bool quantum, bool hessianRequired,
                 std::string structureFile, YAML::Node& yamlNode, std::string yamlSettingsPath, Core::Log& log) {
  // For the objects that we want to be mostly silent
  Core::Log warningLog = Core::Log::silent();
  warningLog.warning.add("cerr", Core::Log::cerrSink());
  warningLog.error.add("cerr", Core::Log::cerrSink());

  if (mode == "calculate" && !quantum) {
    auto mmModel = SwooseUtilities::getChosenMMCalculatorOption(yamlNode);
    auto calculator = getMmCalculatorFromModel(mmModel, manager);
    calculator->setLog(log);
    Utils::nodeToSettings(calculator->settings(), yamlNode, false);
    Utils::PropertyList properties = Utils::Property::Energy | Utils::Property::Gradients;
    if (hessianRequired)
      properties = Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian;
    Tasks::runMMCalculationTask(*calculator, structureFile, properties, log);
  }
  else if (mode == "calculate") {
    std::shared_ptr<Qmmm::QmmmCalculator> calculator;
    SwooseUtilities::fillQmmmCalculatorWithUnderlyingCalculators(calculator, yamlNode);
    calculator->setLog(log);
    Utils::PropertyList properties = Utils::Property::Energy | Utils::Property::Gradients;
    Tasks::runQmmmCalculationTask(*calculator, structureFile, properties, log, yamlNode);
  }
  else if (mode == "parametrize") {
    std::shared_ptr<Core::MMParametrizer> parametrizer;
    try {
      parametrizer = manager.get<Core::MMParametrizer>("SFAM_parametrizer");
    }
    catch (const std::runtime_error& e) {
      throw std::runtime_error("The SFAM parametrizer could not be loaded via the module system.");
    }
    parametrizer->setLog(log);
    Utils::nodeToSettings(parametrizer->settings(), yamlNode, true);
    parametrizer->settings().modifyString(SwooseUtilities::SettingsNames::yamlSettingsFilePath, yamlSettingsPath);
    Tasks::runSFAMParametrizationTask(*parametrizer, structureFile, log);
  }
  else if (mode == "md" && !quantum) {
    auto mmModel = SwooseUtilities::getChosenMMCalculatorOption(yamlNode);
    auto calculator = getMmCalculatorFromModel(mmModel, manager);
    calculator->setLog(warningLog);
    Utils::nodeToSettings(calculator->settings(), yamlNode, true);
    Utils::MolecularDynamics molecularDynamics(*calculator);
    Utils::nodeToSettings(molecularDynamics.settings(), yamlNode, true);
    molecularDynamics.settings().normalizeStringCases(); // convert all option names to lower case letters
    Tasks::runMDSimulationTask(molecularDynamics, structureFile, log);
  }
  else if (mode == "md") {
    std::shared_ptr<Qmmm::QmmmCalculator> calculator;
    calculator->setLog(warningLog);
    SwooseUtilities::fillQmmmCalculatorWithUnderlyingCalculators(calculator, yamlNode);
    Utils::nodeToSettings(calculator->settings(), yamlNode, true);
    Utils::MolecularDynamics molecularDynamics(*calculator);
    Utils::nodeToSettings(molecularDynamics.settings(), yamlNode, true);
    Tasks::runMDSimulationTask(molecularDynamics, structureFile, log);
  }
  else if (mode == "optimize" && !quantum) {
    auto mmModel = SwooseUtilities::getChosenMMCalculatorOption(yamlNode);
    auto calculator = getMmCalculatorFromModel(mmModel, manager);
    calculator->setLog(warningLog);
    Utils::nodeToSettings(calculator->settings(), yamlNode, true);
    auto optimizer = std::make_unique<Utils::GeometryOptimizer<Utils::Bfgs>>(*calculator);
    auto optSettings = optimizer->getSettings();
    Utils::nodeToSettings(optSettings, yamlNode, true);
    optimizer->setSettings(optSettings);
    Tasks::runMMOptimizationTask(*calculator, *optimizer, structureFile, log, yamlNode);
  }
  else if (mode == "optimize") {
    std::shared_ptr<Qmmm::QmmmCalculator> calculator;
    calculator->setLog(warningLog);
    SwooseUtilities::fillQmmmCalculatorWithUnderlyingCalculators(calculator, yamlNode);
    auto optimizer = std::make_unique<Utils::QmmmGeometryOptimizer<Utils::Bfgs>>(*calculator);
    auto optSettings = optimizer->getSettings();
    Utils::nodeToSettings(optSettings, yamlNode, true);
    optimizer->setSettings(optSettings);
    Tasks::runQmmmOptimizationTask<Utils::Bfgs>(*calculator, *optimizer, structureFile, log, yamlNode);
  }
  else if (mode == "select_qm") {
    Tasks::runQmRegionSelectionTask(structureFile, log, yamlNode, yamlSettingsPath);
  }
  else if (mode == "prepare-analyze" || mode == "prepare-protonate" || mode == "prepare-finalize" || mode == "prepare-automate") {
    auto processor = std::make_unique<Scine::StructurePreparation::StructureProcessor>();
    Utils::nodeToSettings(processor->settings(), yamlNode, false); // TODO: allow superfluous?
    Tasks::runPDBPreparationTask(*processor, structureFile, mode, log);
  }
  else {
    throw std::runtime_error(
        "Your specified mode is not valid. Options are: parametrize, calculate, md, optimize, and select_qm.");
  }
}

} // namespace TaskManagement
} // namespace Swoose
} // namespace Scine

#endif // SWOOSE_TASKMANAGEMENT_H
