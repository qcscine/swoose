/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BondedTerm.h"
#include "../MMExceptions.h"
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>

namespace Scine {
namespace MolecularMechanics {

BondedTerm::BondedTerm(AtomIndex firstAtom, AtomIndex secondAtom, const Bond& bond, const BondType& typeOfBond)
  : firstAtom_(firstAtom), secondAtom_(secondAtom), bond_(bond), typeOfBond_(typeOfBond) {
}

BondedTerm::~BondedTerm() = default;

double BondedTerm::evaluateBondTerm(const Utils::PositionCollection& positions,
                                    Utils::AtomicSecondDerivativeCollection& derivatives) const {
  if (this->disabled_)
    return 0.0;
  if (!bond_.hasParameters()) // Check this only if this term is not disabled
    throw MMBondParametersNotAvailableException(typeOfBond_.a1, typeOfBond_.a2);

  auto R = positions.row(secondAtom_) - positions.row(firstAtom_);

  auto v = Utils::AutomaticDifferentiation::get3Dfrom1D<Utils::DerivativeOrder::Two>(bond_.getInteraction(R.norm()), R);
  derivatives[secondAtom_] += v;
  derivatives[firstAtom_] += v.opposite();

  return v.value();
}

BondType BondedTerm::getTypeOfBond() const {
  return typeOfBond_;
}

int BondedTerm::getFirstAtom() const {
  return firstAtom_;
}

int BondedTerm::getSecondAtom() const {
  return secondAtom_;
}

} // namespace MolecularMechanics
} // namespace Scine