/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_HYDROGENBONDEVALUATOR_H
#define MOLECULARMECHANICS_HYDROGENBONDEVALUATOR_H

#include "HydrogenBondTerm.h"

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace MolecularMechanics {
/**
 * @class HydrogenBondEvaluator HydrogenBondEvaluator.h
 * @brief This class evaluates the total energy and its derivative from all hydrogen bonds.
 */
class HydrogenBondEvaluator {
 public:
  /**
   * @brief Constructor from structure and partial atomic charges.
   */
  HydrogenBondEvaluator(const Utils::AtomCollection& structure, const std::vector<double>& atomicCharges);
  /**
   * @brief This function calculates and returns the total hydrogen bond energy and updates the derivatives.
   */
  double evaluate(Utils::AtomicSecondDerivativeCollection& derivatives);
  /**
   * @brief Sets a vector of instances of the HydrogenBondTerm class.
   */
  void setHydrogenBondTerms(std::vector<HydrogenBondTerm>&& hydrogenBondTerms);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::AtomCollection& structure_;
  const std::vector<double>& atomicCharges_;

  std::vector<HydrogenBondTerm> hydrogenBondTerms_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_HYDROGENBONDEVALUATOR_H
