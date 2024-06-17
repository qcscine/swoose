/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ImproperDihedralsEvaluator.h"

namespace Scine {
namespace MolecularMechanics {

ImproperDihedralsEvaluator::ImproperDihedralsEvaluator(const Utils::PositionCollection& positions)
  : positions_(positions) {
}

double ImproperDihedralsEvaluator::evaluate(Utils::AtomicSecondDerivativeCollection& derivatives) {
  double energy = 0.0;

  for (auto& improper : improperDihedrals_) {
    energy += improper.evaluateImproperDihedralTerm(positions_, derivatives);
  }

  return energy;
}

void ImproperDihedralsEvaluator::setImproperDihedralTerms(std::vector<ImproperDihedralTerm>&& improperDihedralTerms) {
  improperDihedrals_ = improperDihedralTerms;
}

} // namespace MolecularMechanics
} // namespace Scine
