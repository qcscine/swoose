/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LennardJones.h"
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>

namespace Scine {
namespace MolecularMechanics {

LennardJones::LennardJones(double aij, double bij) : aij_(aij), bij_(bij), parametersAreAvailable_(true) {
}

LennardJones::LennardJones() : aij_(0.0), bij_(0.0), parametersAreAvailable_(false) {
}

bool LennardJones::hasParameters() const {
  return parametersAreAvailable_;
}

Utils::AutomaticDifferentiation::Second1D LennardJones::getInteraction(double distance) const {
  auto dist2 = distance * distance;
  auto dist3 = dist2 * distance;
  auto dist6 = dist3 * dist3;

  auto dist6inv = Utils::AutomaticDifferentiation::getFromFull<Utils::DerivativeOrder::Two>(
      1.0 / dist6, -6.0 / (dist6 * distance), 42.0 / (dist6 * dist2));
  auto dist12inv = dist6inv * dist6inv;

  return aij_ * dist12inv - bij_ * dist6inv;
}

} // namespace MolecularMechanics
} // namespace Scine