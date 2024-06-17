/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_REFERENCECALCULATIONSIO_H
#define MMPARAMETRIZATION_REFERENCECALCULATIONSIO_H

#include <memory>
#include <vector>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace StructurePreparation {
struct TitrableSite;
} // namespace StructurePreparation

namespace Utils {
class Settings;
namespace ExternalQC {
class TurbomoleMainOutputParser;
} // namespace ExternalQC
} // namespace Utils

namespace MMParametrization {
struct ParametrizationData;
struct TitrationResults;

namespace ReferenceCalculationsIO {

/**
 * @brief This function writes the structures/fragments to disk.
 *        It also writes the information about the charge and spin multiplicity of each structure
 *        as well as (if necessary) information about constrained atoms for the fragment optimization.
 *        It is used during the "write" mode in CalculationManager.
 */
void writeXyzFiles(ParametrizationData& data, std::string referenceDataDir, std::shared_ptr<Utils::Settings> settings);

/**
 * @brief This function reads the reference data from disk.
 *        It is used during the "read" mode in CalculationManager.
 */
void readReferenceDataFromFiles(ParametrizationData& data, TitrationResults& titrationResults,
                                std::string referenceDataDir, std::shared_ptr<Utils::Settings> settings, Core::Log& log);

/**
 * @brief Constructs Turbomole parser and returns it.
 *        This parser will already contain the correct output files directory.
 */
Utils::ExternalQC::TurbomoleMainOutputParser getPreparedTurbomoleParser(const std::string& referenceDataDir, int fragmentIndex);
/**
 * @brief This function writes the structures/fragments for all titrable sites to disk. The additional structures
 * correspond to the respective sites in their non-reference (i. e. charged) protonation state.
 */
void writeAdditionalDataForTitration(ParametrizationData& data, int fragmentIndex, int criticalAtomIndex,
                                     std::string referenceDataDir, std::shared_ptr<Utils::Settings> settings);
/**
 * @brief This function reades the reference data for structures/fragments for all titrable sites from disk.
 */
void saveAdditionalStructuresForTitration(ParametrizationData& data, TitrationResults& results, int fragmentIndex,
                                          std::string referenceDataDir);
/**
 * @brief This function reades the reference data for structures/fragments for all titrable sites from disk.
 */
void parseElectronicEnergiesForTitration(ParametrizationData& data, TitrationResults& results, int fragmentIndex,
                                         std::string referenceDataDir, bool parseTurbomoleOutput,
                                         std::shared_ptr<Utils::Settings> settings);
/**
 * @brief Writes a file with atom indices that should be constrained during the optimization.
 *
 * @param constrainedAtoms The atom indices to be constrained.
 * @param constrainedAtomsFile The filename.
 */
void writeConstrainedAtomsFile(const std::vector<int>& constrainedAtoms, std::string& constrainedAtomsFile);

static constexpr const char* nonRefStateDir = "non_ref_state";
} // namespace ReferenceCalculationsIO
} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_REFERENCECALCULATIONSIO_H
