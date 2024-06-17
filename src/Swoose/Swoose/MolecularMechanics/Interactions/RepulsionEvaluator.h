/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_REPULSIONEVALUATOR_H
#define MOLECULARMECHANICS_REPULSIONEVALUATOR_H

#include "RepulsionTerm.h"
#include <vector>

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace Utils {
class AtomCollection;
class FullSecondDerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {

namespace Defaults {
// Default value optimized as described in the paper: J. Chem. Theory Comput. 2020, 16, 3, 1646-1665.
static constexpr double defaultBetaRepulsionParameter = 7.4;
} // namespace Defaults

/**
 * @class RepulsionEvaluator RepulsionEvaluator.h
 * @brief This class evaluates the Pauli repulsion part of the MM model.
 */
class RepulsionEvaluator {
 public:
  /** @brief Constructor */
  explicit RepulsionEvaluator(const Utils::AtomCollection& structure);
  /**
   * @brief This function evaluates the Pauli repulsion energy and adds its the contribution to the derivatives.
   * @param derivatives The derivatives.
   * @param R0 The matrix of D3 cutoff radii.
   * @return The Pauli repulsion energy.
   */
  double evaluate(Utils::FullSecondDerivativeCollection& derivatives, const Eigen::MatrixXd& R0);
  /**
   * @brief Sets all the repulsion terms.
   */
  void setRepulsionTerms(std::vector<RepulsionTerm>&& repulsionTerms);
  /**
   * @brief Setter for the global Pauli repulsion parameter beta.
   */
  void setBetaRepulsion(const double& betaRepulsion);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::AtomCollection& structure_;
  std::vector<RepulsionTerm> repulsions_;
  double betaRepulsion_ = Defaults::defaultBetaRepulsionParameter;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_REPULSIONEVALUATOR_H
