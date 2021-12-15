/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_BONDSEVALUATOR_H
#define MOLECULARMECHANICS_BONDSEVALUATOR_H

#include "BondedTerm.h"
#include <Utils/Typenames.h>
#include <vector>

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace Utils {
class AtomicSecondDerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {
/**
 * @brief BondsEvaluator BondsEvaluator.h
 * @brief This class evaluates the overall energy and derivatives of bonded interactions.
 */
class BondsEvaluator {
 public:
  /**
   * @brief Constructor from positions.
   */
  explicit BondsEvaluator(const Utils::PositionCollection& positions);
  /**
   * @brief This function evaluates and returns the energy for all bonded interactions and updates the derivatives.
   */
  double evaluate(Utils::AtomicSecondDerivativeCollection& derivatives);
  /**
   * @brief Sets a vector of instances of the BondedTerm class.
   */
  void setBondTerms(std::vector<BondedTerm>&& bondTerms);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::PositionCollection& positions_;
  std::vector<BondedTerm> bonds_;
};

} // namespace MolecularMechanics
} // namespace Scine
#endif // MOLECULARMECHANICS_BONDSEVALUATOR_H