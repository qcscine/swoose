/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_GAFFPARAMETERDEFAULTSPROVIDER_H
#define MOLECULARMECHANICS_GAFFPARAMETERDEFAULTSPROVIDER_H

#include <memory>

namespace Scine {
namespace MolecularMechanics {
class GaffParameters;

/**
 * @class GaffParameterDefaultsProvider GaffParameterDefaultsProvider.h
 * @brief This class processes the default GAFF parameters from the file 'GaffDefaultParameters.h' and returns
 *        the result when the `getParameters()` function is called.
 */
class GaffParameterDefaultsProvider {
 public:
  /**
   * @brief Get the GAFF parameters.
   */
  std::unique_ptr<GaffParameters> getParameters();

 private:
  void processParameters(GaffParameters& parameters);
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_GAFFPARAMETERDEFAULTSPROVIDER_H
