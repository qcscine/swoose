/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_BOND_H
#define MOLECULARMECHANICS_BOND_H

#include <Utils/Math/AutomaticDifferentiation/Second1D.h>

namespace Scine {
namespace MolecularMechanics {
/**
 * @class Bond Bond.h
 * @brief Class treating a bonded interaction, based solely on the bond length. (i.e. in 1 dimension)
 */
class Bond {
 public:
  /**
   * @brief Constructor with equilibrium distance and force constant.
   */
  Bond(double equilibriumDistance, double forceConstant);
  /**
   * @brief Constructor without arguments sets both parameters to zero and records that no parameters
   *        are available for this bond, hence throwing an exception if getInteraction() is called. This
   *        will not throw any error if the corresponding bonded term is disabled (e.g., in a QM/MM calculation).
   */
  Bond();
  /**
   * @brief Returns whether parameters are available for this bond.
   */
  bool hasParameters() const;
  /**
   * @brief Evaluates the energy and its derivatives for one bond.
   */
  Utils::AutomaticDifferentiation::Second1D getInteraction(double bondLength) const;

 private:
  double equilibriumDistance_;
  double forceConstant_;
  bool parametersAreAvailable_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_BOND_H