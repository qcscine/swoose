/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Swoose/MolecularMechanics/InteractionExclusion.h"
#include "Topology/IndexedStructuralTopology.h"

namespace Scine {
namespace MolecularMechanics {

InteractionExclusion::InteractionExclusion(unsigned int nAtoms)
  : exclusions_(Eigen::SparseMatrix<bool>(nAtoms, nAtoms)) {
}

void InteractionExclusion::setExclusions(const Eigen::SparseMatrix<bool>& excludedPairs) {
  if (excludedPairs.cols() != exclusions_.cols() || excludedPairs.rows() != exclusions_.rows()) {
    throw std::runtime_error("Matrix dimensions are incorrect. The number of atoms in the system is"
                             " not allowed to change!");
  }
  exclusions_ = excludedPairs;
}

void InteractionExclusion::addExclusions(const Eigen::SparseMatrix<bool>& excludedPairs) {
  if (excludedPairs.cols() != exclusions_.cols() || excludedPairs.rows() != exclusions_.rows()) {
    throw std::runtime_error("Matrix dimensions are incorrect. The added exclusions must refer to the same"
                             " atom indices as the existing ones!");
  }
  exclusions_ = exclusions_ || excludedPairs;
}

const Eigen::SparseMatrix<bool>& InteractionExclusion::getExclusions() {
  return exclusions_;
}
void InteractionExclusion::addExclusions(const IndexedStructuralTopology& topology) {
  std::vector<Eigen::Triplet<bool>> triplets;
  for (const auto& excluded : topology.getExcludedNonBondedContainer()) {
    if (excluded.atom1 > this->exclusions_.cols() || excluded.atom2 > this->exclusions_.cols()) {
      throw std::runtime_error("Exclusion lists does not fit to the atom indices!");
    }
    triplets.emplace_back(excluded.atom1, excluded.atom2, true);
    triplets.emplace_back(excluded.atom2, excluded.atom1, true);
  }
  Eigen::SparseMatrix<bool> exclusions(this->exclusions_.cols(), this->exclusions_.cols());
  exclusions.setFromTriplets(triplets.begin(), triplets.end());
  this->addExclusions(exclusions);
}
void InteractionExclusion::resetExclusions(unsigned int nAtoms) {
  exclusions_ = Eigen::SparseMatrix<bool>(nAtoms, nAtoms);
}

} // namespace MolecularMechanics
} // namespace Scine
