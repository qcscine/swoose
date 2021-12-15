/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Angle.h"

namespace Scine {
namespace MolecularMechanics {

Angle::Angle(double equilibriumAngle, double forceConstant)
  : equilibriumAngle_(equilibriumAngle), forceConstant_(forceConstant), parametersAreAvailable_(true) {
}

Angle::Angle() : equilibriumAngle_(0.0), forceConstant_(0.0), parametersAreAvailable_(false) {
}

bool Angle::hasParameters() const {
  return parametersAreAvailable_;
}

Utils::AutomaticDifferentiation::Second1D Angle::getInteraction(double angle) const {
  Utils::AutomaticDifferentiation::Second1D thetaDif(angle - equilibriumAngle_, 1, 0);
  return 0.5 * forceConstant_ * thetaDif * thetaDif;
}

} // namespace MolecularMechanics
} // namespace Scine