/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_DIHEDRALSEVALUATOR_H
#define MOLECULARMECHANICS_DIHEDRALSEVALUATOR_H

#include "DihedralTerm.h"
#include <Eigen/Core>
#include <vector>

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace MolecularMechanics {

struct DihedralType;

/**
 * @class DihedralsEvaluator DihedralsEvaluator.h
 * @brief This class calculates the total energy with its derivatives for all dihedral terms of the force field.
 */
class DihedralsEvaluator {
 public:
  /** @brief Constructor from positions. */
  explicit DihedralsEvaluator(const Utils::PositionCollection& positions);
  /** @brief Evaluates the energy and updates the derivatives. */
  double evaluate(Utils::AtomicSecondDerivativeCollection& derivatives);
  /** @brief Sets a vector of the instances of the DihedralTerm class. */
  void setDihedralTerms(std::vector<DihedralTerm>&& dihedralTerms);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::PositionCollection& positions_;
  std::vector<DihedralTerm> dihedrals_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_DIHEDRALSEVALUATOR_H
