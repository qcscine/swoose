/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AngleParameters.h"
#include <Utils/Constants.h>

namespace Scine {
namespace MolecularMechanics {

AngleParameters::AngleParameters(double forceConstant, double equilibriumAngle)
  : forceConstant_(forceConstant), equilibriumAngle_(equilibriumAngle) {
}

Angle AngleParameters::toMMAngle() const {
  return {equilibriumAngle_ * Utils::Constants::rad_per_degree, forceConstant_ * Utils::Constants::hartree_per_kCalPerMol};
}

void AngleParameters::setForceConstant(const double& fc) {
  forceConstant_ = fc;
}

void AngleParameters::setEquilibriumAngle(const double& eqAng) {
  equilibriumAngle_ = eqAng;
}

double AngleParameters::getForceConstant() const {
  return forceConstant_;
}

double AngleParameters::getEquilibriumAngle() const {
  return equilibriumAngle_;
}

} // namespace MolecularMechanics
} // namespace Scine
