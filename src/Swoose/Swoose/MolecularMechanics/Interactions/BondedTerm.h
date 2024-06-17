/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_BONDEDTERM_H
#define MOLECULARMECHANICS_BONDEDTERM_H

#include "../Topology/BondType.h"
#include "Bond.h"
#include "InteractionTermBase.h"
#include <Utils/Typenames.h>
#include <string>

namespace Scine {

namespace Utils {
class AtomicSecondDerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {
/**
 * @class BondedTerm BondedTerm.h
 * @brief Class evaluating harmonic bonded interaction between two atoms.
 */
class BondedTerm : public InteractionTermBase {
 public:
  using AtomIndex = int;
  /** @brief Constructor from two atom indices and instances of Bond and BondType classes. */
  BondedTerm(AtomIndex firstAtom, AtomIndex secondAtom, const Bond& bond, const BondType& typeOfBond);
  /** @brief Destructor. */
  ~BondedTerm();

  /**
   * @brief Evaluates energy contribution and adds the derivatives.
   */
  double evaluateBondTerm(const Utils::PositionCollection& positions, Utils::AtomicSecondDerivativeCollection& derivatives) const;

  /**
   * @brief Getter for the bond type.
   */
  BondType getTypeOfBond() const;

  /**
   * @brief Getter for first atom.
   */
  int getFirstAtom() const;

  /**
   * @brief Getter for second atom.
   */
  int getSecondAtom() const;

 private:
  AtomIndex firstAtom_, secondAtom_;
  Bond bond_;
  BondType typeOfBond_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_BONDEDTERM_H
