/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ImproperDihedral.h"
#include <Utils/Constants.h>
#include <cmath>

namespace Scine {
namespace MolecularMechanics {

// Conversion of the planar group threshold to radians.
static double doubleWellPotentialThreshold = ImproperDihedral::getPlanarGroupThreshold() * Utils::Constants::rad_per_degree;

ImproperDihedral::ImproperDihedral(double forceConstant, double equilibriumAngle)
  : forceConstant_(forceConstant), equilibriumAngle_(equilibriumAngle), parametersAreAvailable_(true) {
}

ImproperDihedral::ImproperDihedral() : forceConstant_(0.0), equilibriumAngle_(0.0), parametersAreAvailable_(false) {
}

bool ImproperDihedral::hasParameters() const {
  return parametersAreAvailable_;
}

// Evaluate the interaction energy and derivatives.
// Distinguish between (near) planar groups and trigonal pyramidal groups.
Utils::AutomaticDifferentiation::Second1D ImproperDihedral::getInteraction(double angle) const {
  if (std::abs(equilibriumAngle_) < doubleWellPotentialThreshold) {
    Utils::AutomaticDifferentiation::Second1D thetaDiff(angle - equilibriumAngle_, 1, 0);
    return forceConstant_ * thetaDiff * thetaDiff;
  }
  else {
    Utils::AutomaticDifferentiation::Second1D theta(angle, 1, 0);
    Utils::AutomaticDifferentiation::Second1D thetaZero(equilibriumAngle_, 0, 0);
    return forceConstant_ * (cos(thetaZero) - cos(theta)) * (cos(thetaZero) - cos(theta));
  }
}

double ImproperDihedral::getPlanarGroupThreshold() {
  return planarGroupThreshold_;
}

} // namespace MolecularMechanics
} // namespace Scine
