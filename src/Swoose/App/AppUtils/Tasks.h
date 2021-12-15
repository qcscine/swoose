/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_TASKS_H
#define SWOOSE_TASKS_H

#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/MMParametrizer.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/GeometryOptimization/QmmmGeometryOptimizer.h>
#include <Utils/IO/Yaml.h>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <ostream>

namespace Scine {

namespace Utils {
class MolecularDynamics;
class XyzStreamHandler;
} // namespace Utils

namespace Swoose {
namespace Tasks {

/// @brief Performs a SFAM or GAFF molecular mechanics calculation.
void runMMCalculationTask(Core::Calculator& calculator, const std::string& structureFile,
                          Utils::PropertyList properties, Core::Log& log);

/// @brief Performs a QM/MM calculation (with SFAM or GAFF as MM).
void runQmmmCalculationTask(Core::Calculator& calculator, const std::string& structureFile,
                            Utils::PropertyList properties, Core::Log& log, const YAML::Node& yamlNode);

/// @brief Performs the parametrization of the SFAM molecular mechanics model.
void runSFAMParametrizationTask(Core::MMParametrizer& parametrizer, const std::string& structureFile, Core::Log& log);

/// @brief Performs an MD simulation.
void runMDSimulationTask(Utils::MolecularDynamics& md, const std::string& structureFile, Core::Log& log);

/// @brief Performs SFAM structure optimization task
void runMMOptimizationTask(Core::Calculator& calculator, Utils::GeometryOptimizerBase& optimizer,
                           const std::string& structureFile, Core::Log& log, const YAML::Node& yamlNode);

/// @brief Performs a QM/MM structure optimization (with SFAM or GAFF as MM).
template<class OptimizerType>
void runQmmmOptimizationTask(Core::Calculator& calculator, Utils::QmmmGeometryOptimizer<OptimizerType>& optimizer,
                             const std::string& structureFile, Core::Log& log, const YAML::Node& yamlNode);

/// @brief Performs the selection of the QM region.
void runQmRegionSelectionTask(const std::string& structureFile, Core::Log& log, YAML::Node& yamlNode,
                              std::string yamlSettingsPath);

} // namespace Tasks
} // namespace Swoose
} // namespace Scine

#endif // SWOOSE_TASKS_H
