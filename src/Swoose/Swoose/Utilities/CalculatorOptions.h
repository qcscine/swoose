/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_CALCULATOROPTIONS_H
#define SWOOSEUTILITIES_CALCULATOROPTIONS_H

#include <Utils/IO/Yaml.h>
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

} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_OPTIONNAMES_H
