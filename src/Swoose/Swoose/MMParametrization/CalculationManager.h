/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_CALCULATIONMANAGER_H
#define MMPARAMETRIZATION_CALCULATIONMANAGER_H

#include <memory>
#include <string>
#include <vector>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace Utils {
class Settings;
} // namespace Utils

namespace MMParametrization {
struct ParametrizationData;

/**
 * @class CalculationManager CalculationManager.h
 * @brief This class handles the reference calculations.
 */
class CalculationManager {
 public:
  /**
   * @brief Constructor.
   */
  CalculationManager(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings, Core::Log& log);
  /**
   * @brief This function calculates the reference data (Hessians and atomic charges for subsystems).
   */
  void calculateReferenceData();

 private:
  // The data used within all MM parametrization classes
  ParametrizationData& data_;
  // The settings
  std::shared_ptr<Utils::Settings> settings_;
  // The logger.
  Core::Log& log_;
  // Directory in which the calculation directories are submitted to
  std::string baseWorkingDirectory_;
  // Mode of the reference data generation (direct, reading, writing, database), default: direct
  std::string mode_;
  // Used Quantum Chemistry software, default: orca
  std::string referenceProgram_;
  // The name of the directory where the calculation files are stored to or loaded from if desired
  std::string referenceDataDirectory_;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_CALCULATIONMANAGER_H
