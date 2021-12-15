/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_ELECTROSTATICEVALUATOR_H
#define MOLECULARMECHANICS_ELECTROSTATICEVALUATOR_H

#include "ElectrostaticTerm.h"
#include <vector>

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace MolecularMechanics {

namespace Defaults {
static constexpr double defaultAtomicChargesScalingFactor = 1.0;
} // namespace Defaults

/**
 * @class ElectrostaticEvaluator ElectrostaticEvaluator.h
 * @brief Class evaluating the total energy and derivatives from the electrostatic interactions.
 */
class ElectrostaticEvaluator {
 public:
  /**
   * @brief Constructor from positions and partial atomic charges.
   */
  ElectrostaticEvaluator(const Utils::PositionCollection& positions, const std::vector<double>& atomicCharges);
  /**
   * @brief Function to evaluate and return the total electrostatic energy and update the derivatives.
   */
  double evaluate(Utils::FullSecondDerivativeCollection& derivatives);
  /**
   * @brief Sets a vector of instances of the ElectrostaticTerm class.
   */
  void setElectrostaticTerms(std::vector<ElectrostaticTerm>&& electrostaticTerms);
  /**
   * @brief Constant accessor for the partial atomic charges.
   */
  const std::vector<double>& getAtomicCharges();

  /**
   * @brief Setter for the scaling factor for each atomic charge.
   */
  void setScalingFactor(const double& scalingFactor);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::PositionCollection& positions_;
  const std::vector<double>& atomicCharges_;
  std::vector<ElectrostaticTerm> electrostaticTerms_;
  // Scaling factor applied to each atomic charge.
  double scalingFactor_ = Defaults::defaultAtomicChargesScalingFactor;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_ELECTROSTATICEVALUATOR_H