/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LennardJonesTerm.h"
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/FullSecondDerivativeCollection.h>

namespace Scine {
namespace MolecularMechanics {

LennardJonesTerm::LennardJonesTerm(AtomIndex firstAtom, AtomIndex secondAtom, const LennardJones& lj,
                                   std::shared_ptr<double> cutoffRadius)
  : firstAtom_(firstAtom), secondAtom_(secondAtom), lj_(lj), cutoffRadius_(cutoffRadius) {
}

double LennardJonesTerm::evaluateLennardJonesTerm(const Utils::PositionCollection& positions,
                                                  Utils::FullSecondDerivativeCollection& derivatives) const {
  if (this->disabled_)
    return 0.0;

  auto R = positions.row(secondAtom_) - positions.row(firstAtom_);

  if (R.norm() > *cutoffRadius_)
    return 0.0;

  auto v = Utils::AutomaticDifferentiation::get3Dfrom1D<Utils::DerivativeOrder::Two>(lj_.getInteraction(R.norm()), R);

  derivatives.addDerivative(firstAtom_, secondAtom_, v);

  return v.value();
}

int LennardJonesTerm::getFirstAtom() const {
  return firstAtom_;
}

int LennardJonesTerm::getSecondAtom() const {
  return secondAtom_;
}

} // namespace MolecularMechanics
} // namespace Scine
