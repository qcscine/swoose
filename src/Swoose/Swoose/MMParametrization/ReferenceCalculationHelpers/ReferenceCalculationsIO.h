/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_REFERENCECALCULATIONSIO_H
#define MMPARAMETRIZATION_REFERENCECALCULATIONSIO_H

#include <memory>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace Utils {
class Settings;
namespace ExternalQC {
class TurbomoleMainOutputParser;
} // namespace ExternalQC
} // namespace Utils

namespace MMParametrization {
struct ParametrizationData;

namespace ReferenceCalculationsIO {

/**
 * @brief This function writes the structures/fragments to disk.
 *        It also writes the information about the charge and spin multiplicity of each structure
 *        as well as (if necessary) information about constrained atoms for the fragment optimization.
 *        It is used during the "write" mode in CalculationManager.
 */
void writeXyzFiles(ParametrizationData& data, std::string referenceDataDir);

/**
 * @brief This function reads the reference data from disk.
 *        It is used during the "read" mode in CalculationManager.
 */
void readReferenceDataFromFiles(ParametrizationData& data, std::string referenceDataDir,
                                std::shared_ptr<Utils::Settings> settings, Core::Log& log);

/**
 * @brief Constructs Turbomole parser and returns it.
 *        This parser will already contain the correct output files directory.
 */
Utils::ExternalQC::TurbomoleMainOutputParser getPreparedTurbomoleParser(const std::string& referenceDataDir, int fragmentIndex);

} // namespace ReferenceCalculationsIO
} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_REFERENCECALCULATIONSIO_H
