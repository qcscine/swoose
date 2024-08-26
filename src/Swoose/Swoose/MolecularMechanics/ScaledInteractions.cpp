/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Swoose/MolecularMechanics/ScaledInteractions.h"
#include "Topology/IndexedStructuralTopology.h"

namespace Scine {
namespace MolecularMechanics {

ScaledInteractions::ScaledInteractions(unsigned int nAtoms)
  : scaledAtomPairs_(Eigen::SparseMatrix<bool>(nAtoms, nAtoms)) {
}

void ScaledInteractions::setScaledInteractionPairs(const Eigen::SparseMatrix<bool>& pairs) {
  if (pairs.cols() != scaledAtomPairs_.cols() || pairs.rows() != scaledAtomPairs_.rows()) {
    throw std::runtime_error("Matrix dimensions are incorrect. The number of atoms in the system is"
                             " not allowed to change!");
  }
  scaledAtomPairs_ = pairs;
}

void ScaledInteractions::addScaledInteractionPairs(const Eigen::SparseMatrix<bool>& pairs) {
  if (pairs.cols() != scaledAtomPairs_.cols() || pairs.rows() != scaledAtomPairs_.rows()) {
    throw std::runtime_error("Matrix dimensions are incorrect. The added parameter scaling must refer to the same"
                             " atom indices as the existing ones!");
  }
  scaledAtomPairs_ = scaledAtomPairs_ || pairs;
}

const Eigen::SparseMatrix<bool>& ScaledInteractions::getScaledInteractionPairs() {
  return scaledAtomPairs_;
}

void ScaledInteractions::setInteractionScalingFactor(const double& factor) {
  scalingFactor_ = factor;
}

double ScaledInteractions::getInteractionScalingFactor() {
  return scalingFactor_;
}

void ScaledInteractions::addScaledInteractionPairs(const IndexedStructuralTopology& topology) {
  std::vector<Eigen::Triplet<bool>> triplets;
  for (const auto& scaled : topology.getScaledNonBondedContainer()) {
    triplets.emplace_back(scaled.atom1, scaled.atom2, true);
    triplets.emplace_back(scaled.atom2, scaled.atom1, true);
  }
  Eigen::SparseMatrix<bool> scaling(this->scaledAtomPairs_.cols(), this->scaledAtomPairs_.cols());
  scaling.setFromTriplets(triplets.begin(), triplets.end());
  this->addScaledInteractionPairs(scaling);
}
void ScaledInteractions::resetScaledInteractions(unsigned int nAtoms) {
  scaledAtomPairs_ = Eigen::SparseMatrix<bool>(nAtoms, nAtoms);
}
} // namespace MolecularMechanics
} // namespace Scine
