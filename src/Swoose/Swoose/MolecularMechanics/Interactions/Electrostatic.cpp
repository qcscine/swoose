/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Electrostatic.h"
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>

namespace Scine {
namespace MolecularMechanics {

Electrostatic::Electrostatic(double scalingFactor) : scalingFactor_(scalingFactor) {
}

Utils::AutomaticDifferentiation::Second1D Electrostatic::getInteraction(double distance, double charge1, double charge2) const {
  auto R = Utils::AutomaticDifferentiation::variableWithUnitDerivative<Utils::DerivativeOrder::Two>(distance);
  auto inverseR = 1.0 / R;

  return inverseR * (charge1 * charge2 * scalingFactor_);
}

} // namespace MolecularMechanics
} // namespace Scine