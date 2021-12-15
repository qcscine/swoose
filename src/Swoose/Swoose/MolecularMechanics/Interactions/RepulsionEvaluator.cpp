/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RepulsionEvaluator.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Math/FullSecondDerivativeCollection.h>

namespace Scine {
namespace MolecularMechanics {

RepulsionEvaluator::RepulsionEvaluator(const Utils::AtomCollection& structure) : structure_(structure) {
}

double RepulsionEvaluator::evaluate(Utils::FullSecondDerivativeCollection& derivatives, const Eigen::MatrixXd& R0) {
  RepulsionParameters repulsionParameters(R0, betaRepulsion_);

  double energy = 0.0;

  for (const auto& repulsion : repulsions_) {
    energy += repulsion.evaluateRepulsionTerm(structure_, derivatives, repulsionParameters);
  }

  return energy;
}

void RepulsionEvaluator::setRepulsionTerms(std::vector<RepulsionTerm>&& repulsionTerms) {
  repulsions_ = repulsionTerms;
}

void RepulsionEvaluator::setBetaRepulsion(const double& betaRepulsion) {
  betaRepulsion_ = betaRepulsion;
}

} // namespace MolecularMechanics
} // namespace Scine