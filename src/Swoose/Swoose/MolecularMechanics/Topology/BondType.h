/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_BONDTYPE_H
#define MOLECULARMECHANICS_BONDTYPE_H

#include <string>
#include <utility>

namespace Scine {
namespace MolecularMechanics {

/**
 * @struct BondType BondType.h
 * @brief  Describes a bond uniquely for given atom types. (useful for maps/unordered_maps)
 *
 *         It is made sure that:
 *            - atoms are ordered correctly.
 *            - two bond types can be compared
 */
struct BondType {
  BondType(std::string atom1, std::string atom2);
  bool operator<(const BondType& rhs) const;
  bool operator==(const BondType& rhs) const;

  std::string a1, a2;
};

inline BondType::BondType(std::string atom1, std::string atom2) : a1(std::move(atom1)), a2(std::move(atom2)) {
  if (a1 > a2)
    std::swap(a1, a2);
}

inline bool BondType::operator<(const BondType& rhs) const {
  if (a1 == rhs.a1)
    return a2 < rhs.a2;
  else
    return a1 < rhs.a1;
}

inline bool BondType::operator==(const BondType& rhs) const {
  if ((a1 == rhs.a1) && (a2 == rhs.a2)) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_BONDTYPE_H