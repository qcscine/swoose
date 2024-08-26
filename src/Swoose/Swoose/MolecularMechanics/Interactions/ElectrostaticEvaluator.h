/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_ELECTROSTATICEVALUATOR_H
#define MOLECULARMECHANICS_ELECTROSTATICEVALUATOR_H

#include "Swoose/MolecularMechanics/InteractionExclusion.h"
#include "Swoose/MolecularMechanics/ScaledInteractions.h"
#include <Utils/Math/DerivativeCollection.h>
#include <Utils/Typenames.h>
#include <memory>
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
class ElectrostaticEvaluator : public InteractionExclusion, public ScaledInteractions {
 public:
  /**
   * @brief Constructor from positions and partial atomic charges.
   */
  ElectrostaticEvaluator(const Utils::PositionCollection& positions, const std::vector<double>& atomicCharges);
  /**
   * @brief Function to evaluate and return the total electrostatic energy and update the derivatives.
   */
  double evaluate(Utils::DerivativeCollection& derivatives);
  /**
   * @brief Constant accessor for the partial atomic charges.
   */
  const std::vector<double>& getAtomicCharges();
  /**
   * @brief Setter for the scaling factor for each atomic charge.
   */
  void setScalingFactor(const double& scalingFactor);
  /**
   * @brief Set the cut off radius for the electrostatic interactions.
   * @param cutOffRadius The cut of radius in atomic units.
   */
  void setCutOffRadius(std::shared_ptr<double> cutOffRadius);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::PositionCollection& positions_;
  const std::vector<double>& atomicCharges_;
  // Scaling factor applied to each atomic charge.
  double scalingFactor_ = Defaults::defaultAtomicChargesScalingFactor;
  /**
   * @brief Calculate the energy and force contribution for a single atom.
   * @param atomIndex    The atom index.
   * @param derivatives  The derivative object to which the gradient contributions are added.
   * @param otherAtoms   Sparse list of atoms for which the interaction is calculated.
   * @param scaling      Interaction scaling.
   * @return The energy contribution.
   */
  double evaluateTermsForAtom(unsigned int atomIndex, Utils::DerivativeCollection& derivatives,
                              const Eigen::SparseVector<bool>& otherAtoms, double scaling);
  ///@brief The cut off radius.
  std::shared_ptr<double> cutOffRadius_ = std::make_shared<double>(std::numeric_limits<double>::infinity());
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_ELECTROSTATICEVALUATOR_H
