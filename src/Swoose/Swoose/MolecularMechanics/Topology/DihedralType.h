/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_DIHEDRALTYPE_H
#define MOLECULARMECHANICS_DIHEDRALTYPE_H

#include <string>
#include <utility>

namespace Scine {
namespace MolecularMechanics {
/**
 * @struct DihedralType DihedralType.h
 * @brief Unique descriptor a dihedral for given atom types. (useful for maps/unordered_maps)
 *
 *        It is made sure that:
 *          - atoms are ordered correctly.
 *          - two dihedral types can be compared
 */
struct DihedralType {
  DihedralType(std::string atom1, std::string atom2, std::string atom3, std::string atom4);
  bool operator<(const DihedralType& rhs) const;
  bool operator==(const DihedralType& rhs) const;

  std::string a1, a2, a3, a4;
};

inline DihedralType::DihedralType(std::string atom1, std::string atom2, std::string atom3, std::string atom4)
  : a1(std::move(atom1)), a2(std::move(atom2)), a3(std::move(atom3)), a4(std::move(atom4)) {
  if (a1 > a4) {
    std::swap(a1, a4);
    std::swap(a2, a3);
  }
  else if (a1 == a4 && a2 > a3) {
    std::swap(a2, a3);
  }
}

inline bool DihedralType::operator<(const DihedralType& rhs) const {
  if (a1 == rhs.a1) {
    if (a2 == rhs.a2) {
      if (a3 == rhs.a3)
        return a4 < rhs.a4;
      else
        return a3 < rhs.a3;
    }
    else
      return a2 < rhs.a2;
  }
  else
    return a1 < rhs.a1;
}

inline bool DihedralType::operator==(const DihedralType& rhs) const {
  if ((a2 == rhs.a2 && a3 == rhs.a3) || (a2 == rhs.a3 && a3 == rhs.a2)) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_DIHEDRALTYPE_H