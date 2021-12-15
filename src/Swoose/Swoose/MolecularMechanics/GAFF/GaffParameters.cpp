/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaffParameters.h"
#include "../AtomTypesHolder.h"
#include "../Interactions/Electrostatic.h"
#include "../MMExceptions.h"
#include <Utils/Geometry/AtomCollection.h>

namespace Scine {
namespace MolecularMechanics {

std::vector<Dihedral> GaffParameters::getMMDihedrals(std::string t1, std::string t2, std::string t3, std::string t4) const {
  std::vector<Dihedral> dihedralsForGivenAtomTypes;

  auto dihedralType = DihedralType(t1, t2, t3, t4);
  bool fullySpecifiedExists = (dihedrals_.find(dihedralType) != dihedrals_.end());

  // If fully specified dihedral doesn't exist, take the one with any end atom types
  if (!fullySpecifiedExists)
    dihedralType = DihedralType("X", t2, t3, "X");

  bool noDihedralTermExists = true;
  for (auto it = dihedrals_.equal_range(dihedralType).first; it != dihedrals_.equal_range(dihedralType).second; ++it) {
    noDihedralTermExists = false;
    // Add only if barrier is non-zero
    if (!(it->second.isZero()))
      dihedralsForGivenAtomTypes.push_back(it->second.toMMDihedral());
  }

  if (noDihedralTermExists)
    throw MMDihedralParametersNotAvailableException(t1, t2, t3, t4);

  return dihedralsForGivenAtomTypes;
}

std::vector<Dihedral> GaffParameters::getMMImproperDihedrals(std::string central, std::string t2, std::string t3,
                                                             std::string t4) const {
  std::vector<Dihedral> improperDihedralsForGivenAtomTypes;

  auto improperDihedralType = ImproperDihedralType(central, t2, t3, t4);

  // If the specified type does not exist, look for others
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, "X", t3, t4);
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, t2, "X", t4);
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, t2, t3, "X");
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, "X", "X", t4);
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, "X", t3, "X");
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, t2, "X", "X");
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, "X", "X", "X");

  for (auto it = improperDihedrals_.equal_range(improperDihedralType).first;
       it != improperDihedrals_.equal_range(improperDihedralType).second; ++it) {
    improperDihedralsForGivenAtomTypes.push_back(it->second.toMMDihedral());
  }

  // NB: no error thrown if empty

  return improperDihedralsForGivenAtomTypes;
}

LennardJones GaffParameters::getMMLennardJones(std::string t1, std::string t2, double scalingFactor) const {
  auto vdwAtom1 = lennardJonesPairs_.find(t1);
  auto vdwAtom2 = lennardJonesPairs_.find(t2);

  if (vdwAtom1 == lennardJonesPairs_.end())
    throw MMLjParametersNotAvailableException(t1);
  if (vdwAtom2 == lennardJonesPairs_.end())
    throw MMLjParametersNotAvailableException(t2);

  return vdwAtom1->second.toMMLennardJones(vdwAtom2->second, scalingFactor);
}

void GaffParameters::addDihedral(DihedralType dihedralType, DihedralParameters dihedralParameters) {
  dihedrals_.emplace(dihedralType, dihedralParameters);
}

void GaffParameters::addImproperDihedral(ImproperDihedralType improperDihedralType,
                                         DihedralParameters improperDihedralParameters) {
  improperDihedrals_.emplace(improperDihedralType, improperDihedralParameters);
}

void GaffParameters::addLennardJones(std::string atomType,
                                     Scine::MolecularMechanics::LennardJonesParameters lennardJonesParameters) {
  lennardJonesPairs_.emplace(atomType, lennardJonesParameters);
}

} // namespace MolecularMechanics
} // namespace Scine
