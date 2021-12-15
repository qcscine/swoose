/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMMMGRADIENTSEVALUATOR_H
#define SWOOSE_QMMM_QMMMGRADIENTSEVALUATOR_H

#include <Utils/Typenames.h>
#include <list>

namespace Scine {

namespace Utils {
class AtomCollection;
} // namespace Utils

namespace Qmmm {

class QmmmGradientsEvaluator {
 public:
  /**
   * @brief Constructor.
   * @param qmGradients Gradients of the QM calculation.
   * @param mmGradients Gradients of the MM calculation.
   * @param pcGradients Gradients Contributions from the point charges.
   *                    This is empty if there is no electrostatic embedding.
   * @param listOfQmAtoms A vector containing the indices of the QM atoms.
   * @param mmBoundaryAtoms A vector that contains the indices of the MM atoms that are bonded to QM atoms.
   * @param listsOfNeighbors A vector containing a list of neighboring atom indices (covalently bonded) for every atom.
   * @param fullStructure The full system's molecular structure.
   * @param qmRegion The QM regions molecular structure containing also the link atoms.
   */
  QmmmGradientsEvaluator(const Utils::GradientCollection& qmGradients, const Utils::GradientCollection& mmGradients,
                         const Utils::GradientCollection& pcGradients, const std::vector<int>& listOfQmAtoms,
                         const std::vector<int>& mmBoundaryAtoms, const std::vector<std::list<int>>& listsOfNeighbors,
                         const Utils::AtomCollection& fullStructure, const Utils::AtomCollection& qmRegion);

  /**
   * @brief Calculates the QM/MM gradient from all the components given in the constructor.
   * @return The QM/MM gradients.
   */
  Utils::GradientCollection calculateQmmmGradients();

 private:
  /*
   * @brief Adds the gradients contributions arising from the QM-MM boundary to the given QM/MM gradients.
   */
  void addBoundaryGradientsContributions(Utils::GradientCollection& qmmmGradients);
  /*
   * @brief Calculates the gradient contribution arising from one specific QM-MM boundary.
   */
  std::pair<double, double>
  calculateGradientContributionForOneBoundary(const Eigen::Ref<Eigen::RowVector3d> qmAtomPosition,
                                              const Eigen::Ref<Eigen::RowVector3d> mmAtomPosition,
                                              const Eigen::Ref<Eigen::RowVector3d> linkAtomPosition,
                                              const Eigen::Ref<Eigen::RowVector3d> linkAtomGradient, int dimension);

  // Several const ref members initiated in the constructor:
  const Utils::GradientCollection& qmGradients_;
  const Utils::GradientCollection& mmGradients_;
  const Utils::GradientCollection& pcGradients_;
  const std::vector<int>& listOfQmAtoms_;
  const std::vector<int>& mmBoundaryAtoms_;
  const std::vector<std::list<int>>& listsOfNeighbors_;
  const Utils::AtomCollection& fullStructure_;
  const Utils::AtomCollection& qmRegion_;
};

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMMMGRADIENTSEVALUATOR_H
