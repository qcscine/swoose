/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ElectrostaticEvaluator.h"

namespace Scine {
namespace MolecularMechanics {
ElectrostaticEvaluator::ElectrostaticEvaluator(const Utils::PositionCollection& positions, const std::vector<double>& atomicCharges)
  : positions_(positions), atomicCharges_(atomicCharges) {
}

double ElectrostaticEvaluator::evaluate(Utils::FullSecondDerivativeCollection& derivatives) {
  double energy = 0.0;

  for (const auto& electrostaticTerm : electrostaticTerms_) {
    energy += electrostaticTerm.evaluateElectrostaticTerm(positions_, derivatives, atomicCharges_, scalingFactor_);
  }

  return energy;
}

void ElectrostaticEvaluator::setElectrostaticTerms(std::vector<ElectrostaticTerm>&& electrostaticTerms) {
  electrostaticTerms_ = electrostaticTerms;
}

const std::vector<double>& ElectrostaticEvaluator::getAtomicCharges() {
  return atomicCharges_;
}

void ElectrostaticEvaluator::setScalingFactor(const double& scalingFactor) {
  scalingFactor_ = scalingFactor;
}

} // namespace MolecularMechanics
} // namespace Scine