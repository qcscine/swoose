/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMIMPROPERDIHEDRALSEVALUATOR_H
#define MMIMPROPERDIHEDRALSEVALUATOR_H

#include "ImproperDihedralTerm.h"
#include <Eigen/Core>
#include <vector>

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace MolecularMechanics {

class ImproperDihedralType;

/**
 * @class ImproperDihedralsEvaluator ImproperDihedralsEvaluator.h
 * @brief This class calculates the total energy and its derivatives for all improper dihedral terms.
 */
class ImproperDihedralsEvaluator {
 public:
  /**
   * @brief Constructor from positions.
   */
  explicit ImproperDihedralsEvaluator(const Utils::PositionCollection& positions);
  /**
   * @brief Evaluates the energy and updates the derivatives.
   */
  double evaluate(Utils::AtomicSecondDerivativeCollection& derivatives);
  /**
   * @brief Sets a vector of instances of the ImproperDihedralTerm class.
   */
  void setImproperDihedralTerms(std::vector<ImproperDihedralTerm>&& improperDihedralTerms);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::PositionCollection& positions_;

  std::vector<ImproperDihedralTerm> improperDihedrals_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MMIMPROPERDIHEDRALSEVALUATOR_H
