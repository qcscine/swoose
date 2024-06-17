/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_PARAMETERFILEWRITER_H
#define MMPARAMETRIZATION_PARAMETERFILEWRITER_H

#include <list>
#include <string>
#include <vector>

namespace Scine {
namespace Utils {
class Settings;
}

namespace MolecularMechanics {
class SfamParameters;
} // namespace MolecularMechanics

namespace MMParametrization {

/**
 * @class ParameterFileWriter ParameterFileWriter.h
 * @brief This class writes the results of a SFAM parametrization (parameters) to file.
 */
class ParameterFileWriter {
 public:
  ParameterFileWriter() = delete;
  /**
   * @brief Writes MM parameters to a file and adds the parametrization settings to the header.
   * @param filename Name of the file to which the MM parameters are written.
   * @param parameters The MM parameters.
   * @param settings The settings.
   */
  static void writeSfamParametersToFile(const std::string& filename, const MolecularMechanics::SfamParameters& parameters,
                                        const Utils::Settings& settings);
  /**
   * @brief Writes MM parameters to a file.
   * @param filename Name of the file to which the MM parameters are written.
   * @param parameters The MM parameters.
   */
  static void writeSfamParametersToFile(const std::string& filename, const MolecularMechanics::SfamParameters& parameters);
  // Implementation of the actual writer.
  static void writeSfamParametersToFileImpl(std::ofstream& parFile, const MolecularMechanics::SfamParameters& parameters);
  // Translates the parametrization settings to a joint string that is written to the header.
  static std::string generateAdditionalInformation(const Utils::Settings& settings);
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_PARAMETERFILEWRITER_H
