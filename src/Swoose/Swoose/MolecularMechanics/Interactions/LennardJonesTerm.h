/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_LENNARDJONESTERM_H
#define MOLECULARMECHANICS_LENNARDJONESTERM_H

#include "InteractionTermBase.h"
#include "LennardJones.h"
#include <Utils/Typenames.h>
#include <memory>
#include <vector>

namespace Scine {

namespace Utils {
class FullSecondDerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {
/**
 * @class LennardJonesTerm LennardJonesTerm.h
 * @brief Class evaluating electrostatic interaction between two atoms.
 */
class LennardJonesTerm : public InteractionTermBase {
 public:
  using AtomIndex = int;

  /**
   * @brief Constructor from two atom indices and an instance of the LennardJones class.
   */
  LennardJonesTerm(AtomIndex firstAtom, AtomIndex secondAtom, const LennardJones& lj, std::shared_ptr<double> cutoffRadius);

  /**
   * @brief Evaluates energy contribution and adds the derivatives.
   */
  double evaluateLennardJonesTerm(const Utils::PositionCollection& positions,
                                  Utils::FullSecondDerivativeCollection& derivatives) const;

  /** @brief Getter for index of first atom. */
  int getFirstAtom() const;
  /** @brief Getter for index of second atom. */
  int getSecondAtom() const;

 private:
  AtomIndex firstAtom_, secondAtom_;
  LennardJones lj_;
  // Pointer so that it is updated whenever the settings are updated.
  // Reference does not work because then the class loses its assignment operator.
  std::shared_ptr<double> cutoffRadius_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_LENNARDJONESTERM_H