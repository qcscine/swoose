/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_HYDROGENBONDTERM_H
#define MOLECULARMECHANICS_HYDROGENBONDTERM_H

#include "HydrogenBond.h"
#include "InteractionTermBase.h"
#include <vector>

namespace Scine {

namespace Utils {
class AtomicSecondDerivativeCollection;
class AtomCollection;
} // namespace Utils

namespace MolecularMechanics {
/**
 * @class HydrogenBondTerm HydrogenBondTerm.h
 * @brief This class evaluates the energy and its derivatives for one hydrogen bond term.
 */
class HydrogenBondTerm : public InteractionTermBase {
 public:
  using AtomIndex = int;
  /**
   * @brief Constructor from the three corresponding atoms and an instance of the HydrogenBond class.
   */
  HydrogenBondTerm(AtomIndex donorAtom, AtomIndex hydrogenAtom, AtomIndex acceptorAtom, const HydrogenBond& hydrogenBond);
  /**
   * @brief Destructor.
   */
  ~HydrogenBondTerm();

  /**
   * @brief Evaluates energy contribution and adds the derivatives.
   */
  double evaluateHydrogenBondTerm(const Utils::AtomCollection& structure, Utils::AtomicSecondDerivativeCollection& derivatives,
                                  const std::vector<double>& atomicCharges) const;

  /// @brief Getter for the donor atom index.
  int getDonorAtom() const;

  /// @brief Getter for the hydrogen atom index.
  int getHydrogenAtom() const;

  /// @brief Getter for the acceptor atom index.
  int getAcceptorAtom() const;

 private:
  AtomIndex donorAtom_, hydrogenAtom_, acceptorAtom_;
  HydrogenBond hydrogenBond_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_HYDROGENBONDTERM_H
