/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ElectrostaticTerm.h"
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/FullSecondDerivativeCollection.h>

namespace Scine {
namespace MolecularMechanics {
ElectrostaticTerm::ElectrostaticTerm(AtomIndex firstAtom, AtomIndex secondAtom, const Electrostatic& electrostatic,
                                     std::shared_ptr<double> cutoffRadius)
  : firstAtom_(firstAtom), secondAtom_(secondAtom), electrostatic_(electrostatic), cutoffRadius_(cutoffRadius) {
}

double ElectrostaticTerm::evaluateElectrostaticTerm(const Utils::PositionCollection& positions,
                                                    Utils::FullSecondDerivativeCollection& derivatives,
                                                    const std::vector<double>& atomicCharges,
                                                    const double& scalingFactorForEachCharge) const {
  if (this->disabled_)
    return 0.0;

  auto R = positions.row(secondAtom_) - positions.row(firstAtom_);

  if (R.norm() > *cutoffRadius_)
    return 0.0;

  auto v = Utils::AutomaticDifferentiation::get3Dfrom1D<Utils::DerivativeOrder::Two>(
      electrostatic_.getInteraction(R.norm(), scalingFactorForEachCharge * atomicCharges[firstAtom_],
                                    scalingFactorForEachCharge * atomicCharges[secondAtom_]),
      R);

  derivatives.addDerivative(firstAtom_, secondAtom_, v);

  return v.value();
}

int ElectrostaticTerm::getFirstAtom() const {
  return firstAtom_;
}

int ElectrostaticTerm::getSecondAtom() const {
  return secondAtom_;
}

} // namespace MolecularMechanics
} // namespace Scine