/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "HydrogenBondEvaluator.h"

namespace Scine {
namespace MolecularMechanics {

HydrogenBondEvaluator::HydrogenBondEvaluator(const Utils::AtomCollection& structure, const std::vector<double>& atomicCharges)
  : structure_(structure), atomicCharges_(atomicCharges) {
}

double HydrogenBondEvaluator::evaluate(Utils::AtomicSecondDerivativeCollection& derivatives) {
  double energy = 0.0;

  for (const auto& hydrogenBondTerm : hydrogenBondTerms_) {
    energy += hydrogenBondTerm.evaluateHydrogenBondTerm(structure_, derivatives, atomicCharges_);
  }

  return energy;
}

void HydrogenBondEvaluator::setHydrogenBondTerms(std::vector<HydrogenBondTerm>&& hydrogenBondTerms) {
  hydrogenBondTerms_ = hydrogenBondTerms;
}

} // namespace MolecularMechanics
} // namespace Scine