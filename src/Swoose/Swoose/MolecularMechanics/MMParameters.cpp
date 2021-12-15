/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MMParameters.h"
#include "MMExceptions.h"

namespace Scine {
namespace MolecularMechanics {

Bond MMParameters::getMMBond(std::string t1, std::string t2) const {
  auto bond_ptr = bonds_.find(BondType(t1, t2));
  if (bond_ptr == bonds_.end())
    throw MMBondParametersNotAvailableException(t1, t2);
  else
    return bond_ptr->second.toMMBond();
}

Angle MMParameters::getMMAngle(std::string t1, std::string t2, std::string t3) const {
  auto angle_ptr = angles_.find(AngleType(t1, t2, t3));
  if (angle_ptr == angles_.end())
    throw MMAngleParametersNotAvailableException(t1, t2, t3);
  else
    return angle_ptr->second.toMMAngle();
}

void MMParameters::addBond(BondType bondType, BondParameters bondParameters) {
  bonds_.emplace(bondType, bondParameters);
}

void MMParameters::addAngle(AngleType angleType, AngleParameters angleParameters) {
  angles_.emplace(angleType, angleParameters);
}

} // namespace MolecularMechanics
} // namespace Scine
