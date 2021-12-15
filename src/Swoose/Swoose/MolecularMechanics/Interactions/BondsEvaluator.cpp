/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BondsEvaluator.h"
#include <Utils/Math/AtomicSecondDerivativeCollection.h>

namespace Scine {
namespace MolecularMechanics {

BondsEvaluator::BondsEvaluator(const Utils::PositionCollection& positions) : positions_(positions) {
}

double BondsEvaluator::evaluate(Utils::AtomicSecondDerivativeCollection& derivatives) {
  double energy = 0.0;

  for (auto& bond : bonds_) {
    energy += bond.evaluateBondTerm(positions_, derivatives);
  }

  return energy;
}

void BondsEvaluator::setBondTerms(std::vector<BondedTerm>&& bondTerms) {
  bonds_ = bondTerms;
}

} // namespace MolecularMechanics
} // namespace Scine