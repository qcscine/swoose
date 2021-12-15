/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_IMPROPERDIHEDRALTYPE_H
#define MOLECULARMECHANICS_IMPROPERDIHEDRALTYPE_H

#include <string>
#include <utility>

namespace Scine {
namespace MolecularMechanics {
/**
 * @struct ImproperDihedralType ImproperDihedralType.h
 * @brief Unique descriptor an improper dihedral for given atom types. (useful for maps/unordered_maps)
 *
 *        It is made sure that:
 *          - atoms are ordered correctly.
 *          - two bond types can be compared
 */
struct ImproperDihedralType {
  ImproperDihedralType(std::string centralAtom, std::string atom2, std::string atom3, std::string atom4);
  bool operator<(const ImproperDihedralType& rhs) const;
  bool operator==(const ImproperDihedralType& rhs) const;

  std::string ac, a2, a3, a4;
};

inline ImproperDihedralType::ImproperDihedralType(std::string centralAtom, std::string atom2, std::string atom3, std::string atom4)
  : ac(std::move(centralAtom)), a2(std::move(atom2)), a3(std::move(atom3)), a4(std::move(atom4)) {
  if (a2 > a4) {
    std::swap(a2, a4);
  }
  if (a2 > a3) {
    std::swap(a2, a3);
  }
  if (a3 > a4) {
    std::swap(a3, a4);
  }
}

inline bool ImproperDihedralType::operator<(const ImproperDihedralType& rhs) const {
  if (ac == rhs.ac) {
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
    return ac < rhs.ac;
}

inline bool ImproperDihedralType::operator==(const ImproperDihedralType& rhs) const {
  if ((ac == rhs.ac) && (a2 == rhs.a2) && (a3 == rhs.a3) && (a4 == rhs.a4)) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_IMPROPERDIHEDRALTYPE_H