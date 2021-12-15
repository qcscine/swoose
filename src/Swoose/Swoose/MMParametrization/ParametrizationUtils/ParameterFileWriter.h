/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_PARAMETERFILEWRITER_H
#define MMPARAMETRIZATION_PARAMETERFILEWRITER_H

#include <list>
#include <string>
#include <vector>

namespace Scine {

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
   * @brief Writes MM parameters to a file.
   * @param filename Name of the file to which the MM parameters are written.
   * @param parameters The MM parameters.
   * @param structure The molecular structure.
   */
  static void writeSfamParametersToFile(const std::string& filename, const MolecularMechanics::SfamParameters& parameters);
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_PARAMETERFILEWRITER_H
