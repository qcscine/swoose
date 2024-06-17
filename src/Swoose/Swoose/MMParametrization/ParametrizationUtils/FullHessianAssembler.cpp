/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FullHessianAssembler.h"
#include "../ParametrizationData.h"
#include <Core/Log.h>
#include <Swoose/Utilities/TopologyUtils.h>

namespace Scine {
namespace MMParametrization {

FullHessianAssembler::FullHessianAssembler(ParametrizationData& data, Core::Log& log) : data_(data), log_(log) {
  fragmentDataDistributor_ = std::make_unique<FragmentDataDistributor>(data, log);
}

void FullHessianAssembler::assembleFullHessian() {
  if (data_.vectorOfStructures.size() == 1) {
    if (data_.vectorOfHessians[0] == nullptr)
      throw std::runtime_error("Fatal error: The reference Hessian could not be obtained.");
    data_.fullHessian = data_.vectorOfHessians[0]->sparseView();
  }
  else {
    data_.fullHessian.resize(3 * data_.numberOfAtoms, 3 * data_.numberOfAtoms);
    assembleFullHessianFromSubsystems();
  }
}

void FullHessianAssembler::assembleFullHessianFromSubsystems() {
  // For bonds
  log_.debug << "Assembling full Hessian - bond components..." << Core::Log::endl;
  assembleBondContributions();
  log_.debug << "Done." << Core::Log::endl;

  // For angles
  log_.debug << "Assembling full Hessian - angle components..." << Core::Log::endl;
  assembleAngleContributions();
  log_.debug << "Done." << Core::Log::endl;

  // For dihedrals
  log_.debug << "Assembling full Hessian - dihedral components..." << Core::Log::endl;
  assembleDihedralContributions();
  log_.debug << "Done." << Core::Log::endl;

  // For improper dihedrals
  log_.debug << "Assembling full Hessian - improper dihedral components..." << Core::Log::endl;
  assembleImproperDihedralContributions();
  log_.debug << "Done." << Core::Log::endl;

  // Print some interesting statistics
  log_.debug << "Density of full Hessian in percent: " << 100.0 * data_.fullHessian.nonZeros() / data_.fullHessian.size()
             << Core::Log::endl;
  double sum = 0.0;
  for (const auto& h : data_.vectorOfHessians) {
    if (h)
      sum += h->size();
  }
  log_.debug << "Number of elements for all subsystem Hessians combined: " << sum << Core::Log::endl;
}

void FullHessianAssembler::assembleBondContributions() {
  for (const auto& b : data_.topology.getBondContainer()) {
    auto relevantAtoms = std::make_pair(b.atom1, b.atom2);
    std::vector<int> orderedCandidatesForHessian = {b.atom1, b.atom2};

    // Add the neighbors of the first candidate atoms to the candidates:
    addNeighboringAtomsToListOfHessianCandidates(orderedCandidatesForHessian);

    addContributionToSparseHessian(relevantAtoms, orderedCandidatesForHessian);
  }
}

void FullHessianAssembler::assembleAngleContributions() {
  for (const auto& a : data_.topology.getAngleContainer()) {
    auto relevantAtoms = std::make_pair(a.atom1, a.atom3);
    std::vector<int> orderedCandidatesForHessian = {a.atom2};

    // Add the neighbors of the first candidate atoms to the candidates:
    addNeighboringAtomsToListOfHessianCandidates(orderedCandidatesForHessian);

    addContributionToSparseHessian(relevantAtoms, orderedCandidatesForHessian);
  }
}

void FullHessianAssembler::assembleDihedralContributions() {
  for (const auto& d : data_.topology.getDihedralContainer()) {
    auto relevantAtoms = std::make_pair(d.atom1, d.atom4);
    std::vector<int> orderedCandidatesForHessian = {d.atom2, d.atom3};

    // Add the neighbors of the first candidate atoms to the candidates:
    addNeighboringAtomsToListOfHessianCandidates(orderedCandidatesForHessian);

    addContributionToSparseHessian(relevantAtoms, orderedCandidatesForHessian);
  }
}

void FullHessianAssembler::assembleImproperDihedralContributions() {
  for (const auto& id : data_.topology.getImproperDihedralContainer()) {
    std::vector<int> orderedCandidatesForHessian = {id.centralAtom, id.atom2, id.atom3, id.atom4};

    // Add the neighbors of the first candidate atoms to the candidates:
    addNeighboringAtomsToListOfHessianCandidates(orderedCandidatesForHessian);

    // For the centralAtom - centralAtom partial Hessian
    auto relevantAtoms1 = std::make_pair(id.centralAtom, id.centralAtom);
    addContributionToSparseHessian(relevantAtoms1, orderedCandidatesForHessian);

    // For the centralAtom - atom2 partial Hessian
    auto relevantAtoms2 = std::make_pair(id.centralAtom, id.atom2);
    addContributionToSparseHessian(relevantAtoms2, orderedCandidatesForHessian);

    // For the centralAtom - atom3 partial Hessian
    auto relevantAtoms3 = std::make_pair(id.centralAtom, id.atom3);
    addContributionToSparseHessian(relevantAtoms3, orderedCandidatesForHessian);

    // For the centralAtom - atom4 partial Hessian
    auto relevantAtoms4 = std::make_pair(id.centralAtom, id.atom4);
    addContributionToSparseHessian(relevantAtoms4, orderedCandidatesForHessian);
  }
}

void FullHessianAssembler::addContributionToSparseHessian(const std::pair<int, int>& relevantAtoms,
                                                          const std::vector<int>& orderedPotentialCandidatesForHessian) {
  for (const auto& candidate : orderedPotentialCandidatesForHessian) {
    if (data_.vectorOfHessians.at(candidate)) {
      const Utils::HessianMatrix& hessian = *data_.vectorOfHessians[candidate];
      auto index1 = std::distance(data_.atomIndexMapping[candidate].begin(),
                                  std::find(data_.atomIndexMapping[candidate].begin(),
                                            data_.atomIndexMapping[candidate].end(), relevantAtoms.first));
      auto index2 = std::distance(data_.atomIndexMapping[candidate].begin(),
                                  std::find(data_.atomIndexMapping[candidate].begin(),
                                            data_.atomIndexMapping[candidate].end(), relevantAtoms.second));
      const Utils::HessianMatrix& subblock = hessian.block<3, 3>(3 * index1, 3 * index2);
      transferSubblockToSparseHessian(subblock, relevantAtoms.first, relevantAtoms.second);
      return;
    }
  }
  throw std::runtime_error("The partial Hessian for the atom pair " + std::to_string(relevantAtoms.first) + " and " +
                           std::to_string(relevantAtoms.second) + " could not be obtained.");
}

// This function is necessary, because there are no block writing operations for sparse matrices in Eigen.
void FullHessianAssembler::transferSubblockToSparseHessian(const Utils::HessianMatrix& subblock, int index1, int index2) {
  assert(subblock.rows() == 3 && subblock.cols() == 3);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      data_.fullHessian.coeffRef(3 * index1 + i, 3 * index2 + j) = subblock(i, j);
    }
  }
}

void FullHessianAssembler::addNeighboringAtomsToListOfHessianCandidates(std::vector<int>& orderedCandidatesForHessian) {
  std::vector<int> firstCandidates = orderedCandidatesForHessian;
  for (const auto& atom : firstCandidates) {
    fragmentDataDistributor_->updateCandidateFragments(atom, orderedCandidatesForHessian);
  }

  // We should have at least 4 candidates to ensure enough redundancy
  if (orderedCandidatesForHessian.size() < 4) {
    for (const auto& atom : firstCandidates) {
      fragmentDataDistributor_->updateCandidateFragmentsWithThirdShellNeighbors(atom, orderedCandidatesForHessian);
    }
  }
}

} // namespace MMParametrization
} // namespace Scine
