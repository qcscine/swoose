/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LennardJonesEvaluator.h"

namespace Scine {
namespace MolecularMechanics {

LennardJonesEvaluator::LennardJonesEvaluator(const Utils::PositionCollection& positions) : positions_(positions) {
}

double LennardJonesEvaluator::evaluate(Utils::FullSecondDerivativeCollection& derivatives) {
  double energy = 0.0;

  for (auto& lj : ljTerms_) {
    energy += lj.evaluateLennardJonesTerm(positions_, derivatives);
  }

  return energy;
}

void LennardJonesEvaluator::setLennardJonesTerms(std::vector<LennardJonesTerm>&& ljTerms) {
  ljTerms_ = ljTerms;
}

} // namespace MolecularMechanics
} // namespace Scine
