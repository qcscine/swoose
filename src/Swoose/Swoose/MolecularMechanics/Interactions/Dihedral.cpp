/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Dihedral.h"

namespace Scine {
namespace MolecularMechanics {

Dihedral::Dihedral(double halfBarrierHeight, int periodicity, double phaseShift)
  : halfBarrierHeight_(halfBarrierHeight), periodicity_(periodicity), phaseShift_(phaseShift), parametersAreAvailable_(true) {
}

Dihedral::Dihedral() : halfBarrierHeight_(0.0), periodicity_(0), phaseShift_(0.0), parametersAreAvailable_(false) {
}

Utils::AutomaticDifferentiation::Second1D Dihedral::getInteraction(double angle) const {
  Utils::AutomaticDifferentiation::Second1D cosTerm(periodicity_ * angle - phaseShift_, periodicity_, 0);
  return halfBarrierHeight_ * (Utils::AutomaticDifferentiation::Second1D(1, 0, 0) + cosPreFactor_ * cos(cosTerm));
}

bool Dihedral::hasParameters() const {
  return parametersAreAvailable_;
}

void Dihedral::setCosinePreFactor(double cosPreFactor) {
  cosPreFactor_ = cosPreFactor;
}

} // namespace MolecularMechanics
} // namespace Scine
