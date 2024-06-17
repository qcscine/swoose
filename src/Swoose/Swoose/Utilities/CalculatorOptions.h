/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_CALCULATOROPTIONS_H
#define SWOOSEUTILITIES_CALCULATOROPTIONS_H

#include <Swoose/QMMM/QmmmCalculator.h>
#include <Utils/CalculatorBasics/CalculationRoutines.h>
#include <Utils/IO/Yaml.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <yaml-cpp/yaml.h>

namespace Scine {
namespace SwooseUtilities {

namespace {
constexpr const char* qmModelKey = "qm_model";
constexpr const char* qmModuleKey = "qm_module";
constexpr const char* mmModelKey = "mm_model";
constexpr const char* defaultMolecularMechanicsModel = "SFAM";
} // namespace

/**
 * @brief Extracts the MM model from the app settings. Used for MM-only and for QM/MM calculations.
 * @param yamlNode The app settings read from the yaml settings file.
 * @return A string that represents the MM model that should be chosen. There is a default, i.e., one does not
 *         have to specify this setting.
 */
inline std::string getChosenMMCalculatorOption(YAML::Node& yamlNode) {
  std::string mmModel = defaultMolecularMechanicsModel;
  if (yamlNode[mmModelKey]) {
    mmModel = yamlNode[mmModelKey].as<std::string>();
    yamlNode.remove(mmModelKey);
  }
  // capitalize MM model string and return it
  std::transform(mmModel.begin(), mmModel.end(), mmModel.begin(), ::toupper);
  return mmModel;
}

/**
 * @brief Extracts the QM model and module for QM/MM calculations from the app settings.
 * @param yamlNode The app settings read from the yaml settings file.
 * @return A pair that contains the QM model (first value) and the SCINE module (second value) which provides it.
 */
inline std::pair<std::string, std::string> getChosenQmCalculatorOption(YAML::Node& yamlNode) {
  std::string qmModel;
  std::string qmModule;
  if (yamlNode[qmModelKey]) {
    qmModel = yamlNode[qmModelKey].as<std::string>();
    yamlNode.remove(qmModelKey);
  }
  if (yamlNode[qmModuleKey]) {
    qmModule = yamlNode[qmModuleKey].as<std::string>();
    yamlNode.remove(qmModuleKey);
  }
  if (qmModel.empty() || qmModule.empty())
    throw std::runtime_error("For a QM/MM calculation, one needs to specify the QM model and module.");

  // Capitalize the QM model:
  std::transform(qmModel.begin(), qmModel.end(), qmModel.begin(), ::toupper);

  // Convert module to only first letter as uppercase:
  std::transform(qmModule.begin(), qmModule.end(), qmModule.begin(), ::tolower);
  qmModule.at(0) = std::toupper(qmModule.at(0));

  return {qmModel, qmModule};
}

/**
 * @brief Helper function that constructs the correct MM and QM calculators from the YAML settings and passes them
 *        to the provided QM/MM calculator, which is passed by reference.
 * @param qmmmCalculator The QM/MM calculator.
 * @param manager The module manager.
 * @param yamlNode The settings as a YAML node.
 */
inline void fillQmmmCalculatorWithUnderlyingCalculators(std::shared_ptr<Scine::Qmmm::QmmmCalculator>& qmmmCalculator,
                                                        YAML::Node& yamlNode) {
  std::string methodFamilies, programs;
  // legacy support, remove in future and only use method_family and program and load QM/MM calc directly
  if (yamlNode[qmModelKey]) {
    methodFamilies = yamlNode[qmModelKey].as<std::string>() + '/' + defaultMolecularMechanicsModel;
  }
  else if (yamlNode[Utils::SettingsNames::methodFamily]) {
    methodFamilies = yamlNode[Utils::SettingsNames::methodFamily].as<std::string>();
    if (methodFamilies.find('/') == std::string::npos) {
      methodFamilies += "/" + std::string(defaultMolecularMechanicsModel);
    }
  }
  else {
    throw std::runtime_error("Please provide two method families for a QM/MM calculation.\n"
                             "For instance, method_family: PM6/SFAM and program: Sparrow/Swoose");
  }
  if (yamlNode[qmModuleKey]) {
    programs = yamlNode[qmModuleKey].as<std::string>() + '/' + "Swoose";
  }
  else if (yamlNode[Utils::SettingsNames::program]) {
    programs = yamlNode[Utils::SettingsNames::program].as<std::string>();
    if (programs.find('/') == std::string::npos) {
      programs += "/" + std::string("Swoose");
    }
  }
  else {
    programs = "Any/Swoose";
  }
  std::string error =
      "The QM/MM calculator could not be loaded via the module system.\nCheck: (i) whether the requested "
      "calculators are available (see manual), (ii) that you have installed all "
      "relevant modules, and (iii) that all necessary environment variables are set (e.g., ORCA_BINARY_PATH or "
      "TURBODIR).";
  try {
    auto loadedCalc = Utils::CalculationRoutines::getCalculator(methodFamilies, programs);
    if (!loadedCalc)
      throw std::runtime_error(error);
    auto castedCalc = std::dynamic_pointer_cast<Qmmm::QmmmCalculator>(loadedCalc);
    if (!castedCalc)
      throw std::runtime_error("The loaded calculator is not a QM/MM calculator.");
    qmmmCalculator = castedCalc;
  }
  catch (const std::runtime_error& e) {
    std::cerr << e.what() << std::endl;
    throw std::runtime_error(error);
  }
}

} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_CALCULATOROPTIONS_H
