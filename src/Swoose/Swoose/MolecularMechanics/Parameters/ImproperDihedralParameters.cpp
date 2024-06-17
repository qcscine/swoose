/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ImproperDihedralParameters.h"
#include <Utils/Constants.h>

namespace Scine {
namespace MolecularMechanics {

ImproperDihedralParameters::ImproperDihedralParameters(double forceConstant, double equilibriumAngle)
  : forceConstant_(forceConstant), equilibriumAngle_(equilibriumAngle) {
}

ImproperDihedral ImproperDihedralParameters::toMMImproperDihedral() const {
  return {forceConstant_ * Utils::Constants::hartree_per_kCalPerMol, equilibriumAngle_ * Utils::Constants::rad_per_degree};
}

void ImproperDihedralParameters::setForceConstant(const double& fc) {
  forceConstant_ = fc;
}

void ImproperDihedralParameters::setEquilibriumAngle(const double& eqAng) {
  equilibriumAngle_ = eqAng;
}

double ImproperDihedralParameters::getForceConstant() const {
  return forceConstant_;
}

double ImproperDihedralParameters::getEquilibriumAngle() const {
  return equilibriumAngle_;
}

} // namespace MolecularMechanics
} // namespace Scine
