/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DihedralsEvaluator.h"
#include <iostream>

namespace Scine {
namespace MolecularMechanics {

DihedralsEvaluator::DihedralsEvaluator(const Utils::PositionCollection& positions) : positions_(positions) {
}

double DihedralsEvaluator::evaluate(Utils::AtomicSecondDerivativeCollection& derivatives) {
  double energy = 0.0;

  for (auto& dihedral : dihedrals_) {
    energy += dihedral.evaluateDihedralTerm(positions_, derivatives);
  }

  return energy;
}

void DihedralsEvaluator::setDihedralTerms(std::vector<DihedralTerm>&& dihedralTerms) {
  dihedrals_ = dihedralTerms;
}

} // namespace MolecularMechanics
} // namespace Scine