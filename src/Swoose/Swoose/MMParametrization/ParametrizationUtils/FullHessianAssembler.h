/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_FULLHESSIANASSEMBLER_H
#define MMPARAMETRIZATION_FULLHESSIANASSEMBLER_H

#include "FragmentDataDistributor.h"
#include <Utils/Typenames.h>
#include <memory>
#include <vector>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace MMParametrization {
struct ParametrizationData;

/**
 * @class FullHessianAssembler FullHessianAssembler.h
 * @brief Assembles the full system's Hessian matrix from the subsystem's Hessians.
 */
class FullHessianAssembler {
 public:
  /**
   * @brief Constructor.
   */
  FullHessianAssembler(ParametrizationData& data, Core::Log& log);
  /**
   * @brief This function takes all the Hessians from the subsystems and the mapping of the parameters
   *        to the subsystems to generate the sparse matrix which is the Hessian of the full system
   *        in the next steps.
   */
  void assembleFullHessian();

 private:
  // Implementation of the public function "assembleFullHessian" for the case of more than one subsystem
  void assembleFullHessianFromSubsystems();
  // Assembles the relevant subblocks of the Hessian for the bonds
  void assembleBondContributions();
  // Assembles the relevant subblocks of the Hessian for the angles
  void assembleAngleContributions();
  // Assembles the relevant subblocks of the Hessian for the dihedals
  void assembleDihedralContributions();
  // Assembles the relevant subblocks of the Hessian for the improper dihedrals
  void assembleImproperDihedralContributions();
  /*
   * This function finds the correct submatrix for a given MM parameter (by relevant atoms)
   * from a given set of candidate subsystems and transfers it to the sparse full Hessian matrix
   */
  void addContributionToSparseHessian(const std::pair<int, int>& relevantAtoms,
                                      const std::vector<int>& orderedPotentialCandidatesForHessian);
  // Transfers a 3x3 matrix to its correct position in the full sparse Hessian matrix
  void transferSubblockToSparseHessian(const Utils::HessianMatrix& subblock, int index1, int index2);
  /*
   * Adds the neighbors of the atoms that are already part of a given MM parameter
   * to the list of candidates for Hessians, too.
   */
  void addNeighboringAtomsToListOfHessianCandidates(std::vector<int>& orderedCandidatesForHessian);
  // The data used within all MM parametrization classes
  ParametrizationData& data_;
  // The logger.
  Core::Log& log_;
  // Pointer to an instance of the fragment data distributor class needed to find fragment candidates to get data from.
  std::unique_ptr<FragmentDataDistributor> fragmentDataDistributor_;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_FULLHESSIANASSEMBLER_H
