/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_DIRECTCALCULATIONSHELPER_H
#define MMPARAMETRIZATION_DIRECTCALCULATIONSHELPER_H

#include <memory>

namespace Scine {

namespace Core {
class Log;
} // namespace Core

namespace Utils {
class Settings;
} // namespace Utils

namespace MMParametrization {
struct ParametrizationData;

namespace DirectCalculationsHelper {

/**
 * @brief Optimizes the molecular structure, calculates the reference Hessian and Mayer bond orders
 *        directly through the SCINE interface.
 *        If the use_gaussian keyword is set to false, it also calculates the atomic partial charges.
 *        It is used during the "direct" mode in CalculationManager.
 */
void performReferenceCalculations(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings,
                                  std::string baseWorkingDir, Core::Log& log);

/**
 * @brief Calculate the atomic charges directly through the SCINE ExternalQC/Gaussian interface
 *        if setting use_gaussian is set to true.
 *        It is used during the "direct" mode in CalculationManager.
 */
void calculateAtomicChargesWithGaussian(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings,
                                        std::string baseWorkingDir, Core::Log& log);

} // namespace DirectCalculationsHelper
} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_DIRECTCALCULATIONSHELPER_H
