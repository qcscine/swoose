/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RepulsionTerm.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/FullSecondDerivativeCollection.h>

namespace Scine {
namespace MolecularMechanics {

RepulsionTerm::RepulsionTerm(AtomIndex firstAtom, AtomIndex secondAtom, const Repulsion& repulsion,
                             std::shared_ptr<double> cutoffRadius)
  : firstAtom_(firstAtom), secondAtom_(secondAtom), repulsion_(repulsion), cutoffRadius_(cutoffRadius) {
}

RepulsionTerm::~RepulsionTerm() = default;

double RepulsionTerm::evaluateRepulsionTerm(const Utils::AtomCollection& structure,
                                            Utils::FullSecondDerivativeCollection& derivatives,
                                            RepulsionParameters& repulsionParameters) const {
  if (this->disabled_)
    return 0.0;

  const Utils::PositionCollection& positions = structure.getPositions();
  const auto& positionA = positions.row(firstAtom_);
  const auto& positionB = positions.row(secondAtom_);
  auto R = positionB - positionA;

  if (R.norm() > *cutoffRadius_)
    return 0.0;

  const auto& elements = structure.getElements();
  const auto& elementA = elements[firstAtom_];
  const auto& elementB = elements[secondAtom_];

  const auto& effectiveChargeA = repulsionParameters.getEffectiveCharge(elementA);
  const auto& effectiveChargeB = repulsionParameters.getEffectiveCharge(elementB);
  const auto& betaRepulsion = repulsionParameters.getBetaRepulsion();

  auto R0 = repulsionParameters.getR0(firstAtom_, secondAtom_);

  auto v = Utils::AutomaticDifferentiation::get3Dfrom1D<Utils::DerivativeOrder::Two>(
      repulsion_.getInteraction(R.norm(), effectiveChargeA, effectiveChargeB, betaRepulsion, R0), R);

  derivatives.addDerivative(firstAtom_, secondAtom_, v);

  return v.value();
}

int RepulsionTerm::getFirstAtom() const {
  return firstAtom_;
}

int RepulsionTerm::getSecondAtom() const {
  return secondAtom_;
}

} // namespace MolecularMechanics
} // namespace Scine