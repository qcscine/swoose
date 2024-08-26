/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_LENNARDJONESEVALUATOR_H
#define MOLECULARMECHANICS_LENNARDJONESEVALUATOR_H

#include "Swoose/MolecularMechanics/InteractionExclusion.h"
#include "Swoose/MolecularMechanics/ScaledInteractions.h"
#include <Utils/Typenames.h>
#include <memory>
#include <vector>

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace Utils {
class DerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {
class GaffParameters;
class AtomTypesHolder;
/**
 * @brief LennardJonesEvaluator LennardJonesEvaluator.h
 * @brief This class evaluates the overall energy and derivatives of Lennard-Jones interactions.
 */
class LennardJonesEvaluator : public InteractionExclusion, public ScaledInteractions {
 public:
  /**
   * @brief Constructor from positions.
   */
  explicit LennardJonesEvaluator(const Utils::PositionCollection& positions);
  /**
   * @brief This function evaluates and returns the energy for all LJ interactions and updates the derivatives.
   */
  double evaluate(Utils::DerivativeCollection& derivatives);
  /**
   * @brief Set the cut off radius for the electrostatic interactions.
   * @param cutOffRadius The cut of radius in atomic units.
   */
  void setCutOffRadius(std::shared_ptr<double> cutOffRadius);
  /**
   * @brief Setter for the parameter object.
   * @param parameters The parameters.
   */
  void setParameters(std::shared_ptr<GaffParameters> parameters);
  /**
   * @brief Setter for the atom types.
   * @param atomTypesHolder The atom tyoes.
   */
  void setAtomTypesHolder(std::shared_ptr<AtomTypesHolder> atomTypesHolder);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::PositionCollection& positions_;

  ///@brief The Parameters.
  std::shared_ptr<GaffParameters> parameters_ = nullptr;
  ///@brief The atom types.
  std::shared_ptr<AtomTypesHolder> atomTypesHolder_ = nullptr;
  ///@brief The cut off radius.
  std::shared_ptr<double> cutOffRadius_ = std::make_shared<double>(std::numeric_limits<double>::infinity());
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
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_LENNARDJONESEVALUATOR_H
