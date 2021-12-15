/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_ANGLESEVALUATOR_H
#define MOLECULARMECHANICS_ANGLESEVALUATOR_H

#include "AngleTerm.h"
#include <vector>

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace MolecularMechanics {

class AngleType;

class AnglesEvaluator {
 public:
  /**
   * @brief Constructor from positions.
   */
  explicit AnglesEvaluator(const Utils::PositionCollection& positions);
  /**
   * @brief This function evaluates and returns the energy for all angle interactions and updates the derivatives.
   */
  double evaluate(Utils::AtomicSecondDerivativeCollection& derivatives);
  /**
   * @brief Sets a vector of instances of AngleTerm.
   */
  void setAngleTerms(std::vector<AngleTerm>&& angleTerms);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::PositionCollection& positions_;
  std::vector<AngleTerm> angles_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_ANGLESEVALUATOR_H