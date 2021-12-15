/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_BASICJOBSUBMISSIONHELPER_H
#define MMPARAMETRIZATION_BASICJOBSUBMISSIONHELPER_H

#include <string>

namespace Scine {

namespace Utils {
class Settings;
} // namespace Utils

namespace MMParametrization {
namespace BasicJobSubmissionHelper {

/**
 * @brief Determines the method family for the employed methods/programs.
 */
std::string determineMethodFamily(std::string method, std::string referenceProgram);
/**
 * @brief Determines whether the basis set for the employed methods/programs should be empty. If not, it returns
 *        the basis set from the settings (given as an argument).
 */
std::string determineBasisSet(const Utils::Settings& settings);

} // namespace BasicJobSubmissionHelper
} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_BASICJOBSUBMISSIONHELPER_H
