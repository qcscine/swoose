/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_LENNARDJONESEVALUATOR_H
#define MOLECULARMECHANICS_LENNARDJONESEVALUATOR_H

#include "LennardJonesTerm.h"
#include <Utils/Typenames.h>
#include <vector>

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace Utils {
class FullSecondDerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {
/**
 * @brief LennardJonesEvaluator LennardJonesEvaluator.h
 * @brief This class evaluates the overall energy and derivatives of Lennard-Jones interactions.
 */
class LennardJonesEvaluator {
 public:
  /**
   * @brief Constructor from positions.
   */
  explicit LennardJonesEvaluator(const Utils::PositionCollection& positions);
  /**
   * @brief This function evaluates and returns the energy for all LJ interactions and updates the derivatives.
   */
  double evaluate(Utils::FullSecondDerivativeCollection& derivatives);
  /**
   * @brief Sets a vector of instances of the LennardJonesTerm class.
   */
  void setLennardJonesTerms(std::vector<LennardJonesTerm>&& ljTerms);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::PositionCollection& positions_;
  std::vector<LennardJonesTerm> ljTerms_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_LENNARDJONESEVALUATOR_H