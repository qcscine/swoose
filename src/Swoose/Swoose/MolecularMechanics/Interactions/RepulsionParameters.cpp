/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RepulsionParameters.h"

namespace Scine {
namespace MolecularMechanics {

RepulsionParameters::RepulsionParameters(const Eigen::MatrixXd& R0, const double& betaRepulsion)
  : R0_(R0), betaRepulsion_(betaRepulsion) {
}

double RepulsionParameters::getR0(int atom1Index, int atom2Index) const {
  return R0_(atom1Index, atom2Index);
}

double RepulsionParameters::getBetaRepulsion() {
  return betaRepulsion_;
}

double RepulsionParameterHelper::getEffectiveCharge(Utils::ElementType element) {
  auto elementZ = Utils::ElementInfo::Z(element);
  int valenceElectronNumbers = valenceElectrons_[elementZ - 1];
  double valenceElectronScalingFactor = valenceElectronScalingFactors_[elementZ - 1];
  return valenceElectronNumbers * valenceElectronScalingFactor;
}

} // namespace MolecularMechanics
} // namespace Scine
