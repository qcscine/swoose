/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RepulsionParameters.h"
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace MolecularMechanics {

RepulsionParameters::RepulsionParameters(const Eigen::MatrixXd& R0, const double& betaRepulsion)
  : R0_(R0), betaRepulsion_(betaRepulsion) {
}

double RepulsionParameters::getEffectiveCharge(Utils::ElementType element) const {
  return getValenceElectronNumbers(element) * getValenceElectronScalingFactor(element);
}

double RepulsionParameters::getR0(int atom1Index, int atom2Index) const {
  return R0_(atom1Index, atom2Index);
}

double RepulsionParameters::getBetaRepulsion() {
  return betaRepulsion_;
}

int RepulsionParameters::getValenceElectronNumbers(Utils::ElementType element) const {
  auto elementZ = Utils::ElementInfo::Z(element);
  return valenceElectrons_[elementZ - 1];
}

double RepulsionParameters::getValenceElectronScalingFactor(Utils::ElementType element) const {
  auto elementZ = Utils::ElementInfo::Z(element);
  return valenceElectronScalingFactors_[elementZ - 1];
}

} // namespace MolecularMechanics
} // namespace Scine
