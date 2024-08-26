/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_INTERACTIONEXCLUSION_H
#define SWOOSE_INTERACTIONEXCLUSION_H

#include <Eigen/Sparse>

namespace Scine {
namespace MolecularMechanics {
class IndexedStructuralTopology;

/**
 * @class InteractionExclusion InteractionExclusion.h
 * @brief This class provides a framework to store exclusions of interactions between atoms.
 *        Such exclusions are required for the non-covalent interactions between bonded atoms or
 *        to eliminate contributions to the energy for the interaction between QM atoms in QM/MM.
 *
 *        The exclusions are stored as a sparse matrix atom x atom of booleans. If the entry is true
 *        for two atoms, their interaction is excluded.
 */
class InteractionExclusion {
 public:
  /**
   * @brief Constructor.
   * @param nAtoms The number of atoms.
   *
   * Note that the number of atoms may be changed later by calling resetExclusions(...).
   */
  InteractionExclusion(unsigned int nAtoms);
  /**
   * @brief Setter for the exclusion matrix.
   * @param excludedPairs The exclusion matrix.
   */
  void setExclusions(const Eigen::SparseMatrix<bool>& excludedPairs);
  /**
   * @brief Add new exclusions.
   * @param excludedPairs The exclusions.
   */
  void addExclusions(const Eigen::SparseMatrix<bool>& excludedPairs);
  /**
   * @brief Add the exclusions encoded in the topology, e.g., 1-3 and 1-2 interactions.
   * @param topology The topology.
   */
  void addExclusions(const IndexedStructuralTopology& topology);
  /**
   * @brief Getter for the exclusions.
   * @return The exclusions.
   */
  const Eigen::SparseMatrix<bool>& getExclusions();
  /**
   * @brief Resize the exclusion matrix to the given number of atoms and remove all existing exlcuions.
   * @param nAtoms The new number of atoms.
   */
  void resetExclusions(unsigned int nAtoms);

 private:
  Eigen::SparseMatrix<bool> exclusions_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // SWOOSE_INTERACTIONEXCLUSION_H
