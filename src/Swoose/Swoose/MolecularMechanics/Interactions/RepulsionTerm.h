/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_REPULSIONTERM_H
#define MOLECULARMECHANICS_REPULSIONTERM_H

#include "InteractionTermBase.h"
#include "Repulsion.h"
#include "RepulsionParameters.h"

namespace Scine {

namespace Utils {
class AtomCollection;
class FullSecondDerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {

/**
 * @class RepulsionTerm RepulsionTerm.h
 * @brief Class evaluating repulsive non-bonded interaction between two atoms.
 */
class RepulsionTerm : public InteractionTermBase {
 public:
  using AtomIndex = int;

  /**
   * @brief Constructor
   * @param firstAtom Index of first atom.
   * @param secondAtom Index of second atom.
   * @param repulsion The corresponding instance of the Repulsion class.
   */
  RepulsionTerm(AtomIndex firstAtom, AtomIndex secondAtom, const Repulsion& repulsion, std::shared_ptr<double> cutoffRadius);
  /** @brief Destructor */
  ~RepulsionTerm();

  /**
   * @brief Evaluates the energy contribution and adds the derivatives.
   */
  double evaluateRepulsionTerm(const Utils::AtomCollection& structure, Utils::FullSecondDerivativeCollection& derivatives,
                               RepulsionParameters& repulsionParameters) const;

  /** @brief Getter for index of first atom. */
  int getFirstAtom() const;
  /** @brief Getter for index of second atom. */
  int getSecondAtom() const;

 private:
  AtomIndex firstAtom_, secondAtom_;
  Repulsion repulsion_;
  // Pointer so that it is updated whenever the settings are updated.
  // Reference does not work because then the class loses its assignment operator.
  std::shared_ptr<double> cutoffRadius_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_REPULSIONTERM_H