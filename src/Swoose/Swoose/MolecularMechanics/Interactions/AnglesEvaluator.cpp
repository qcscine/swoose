/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AnglesEvaluator.h"

namespace Scine {
namespace MolecularMechanics {
AnglesEvaluator::AnglesEvaluator(const Utils::PositionCollection& positions) : positions_(positions) {
}

double AnglesEvaluator::evaluate(Utils::AtomicSecondDerivativeCollection& derivatives) {
  double energy = 0.0;

  for (auto& angle : angles_) {
    energy += angle.evaluateAngleTerm(positions_, derivatives);
  }

  return energy;
}

void AnglesEvaluator::setAngleTerms(std::vector<AngleTerm>&& angleTerms) {
  angles_ = angleTerms;
}

} // namespace MolecularMechanics
} // namespace Scine
