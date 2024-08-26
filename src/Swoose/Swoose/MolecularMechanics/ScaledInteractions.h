/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_SCALEDINTERACTIONS_H
#define SWOOSE_SCALEDINTERACTIONS_H

#include <Eigen/Sparse>

namespace Scine {
namespace MolecularMechanics {
class IndexedStructuralTopology;

/**
 * @class ScaledInteractions ScaledInteractions.h
 * @brief This class keeps track which pair-wise interaction between atoms is scaled.
 *        The scaling information is stored as a sparse matrix atoms x atoms of booleans.
 */
class ScaledInteractions {
 public:
  /**
   * @brief Constructor.
   * @param nAtoms The number of atoms.
   */
  ScaledInteractions(unsigned int nAtoms);
  /**
   * @brief Setter for the scaling information.
   * @param pairs The scaled atom pairs.
   */
  void setScaledInteractionPairs(const Eigen::SparseMatrix<bool>& pairs);
  /**
   * @brief Add new scaled interactions.
   * @param pairs The additional scaled atom pairs.
   */
  void addScaledInteractionPairs(const Eigen::SparseMatrix<bool>& pairs);
  /**
   * @brief Add new scaled interactions from the topology, e.g., 1-4 interactions.
   * @param topology The topology.
   */
  void addScaledInteractionPairs(const IndexedStructuralTopology& topology);
  /**
   * @brief Getter for scaling information.
   * @return The scaling matrix.
   */
  const Eigen::SparseMatrix<bool>& getScaledInteractionPairs();
  /**
   * @brief Setter for the interaction scaling factor.
   * @param factor The factor.
   */
  void setInteractionScalingFactor(const double& factor);
  /**
   * @brief Getter for the interaction scaling factor.
   * @return The interaction scaling factor.
   */
  double getInteractionScalingFactor();
  /**
   * @brief Resize the scaling information matrix to the nem number of atoms.
   * @param nAtoms The number of atoms.
   */
  void resetScaledInteractions(unsigned int nAtoms);

 private:
  Eigen::SparseMatrix<bool> scaledAtomPairs_;
  double scalingFactor_ = 0.5;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // SWOOSE_SCALEDINTERACTIONS_H
