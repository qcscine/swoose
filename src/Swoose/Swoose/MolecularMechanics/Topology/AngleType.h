/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_ANGLETYPE_H
#define MOLECULARMECHANICS_ANGLETYPE_H

#include <string>
#include <utility>

namespace Scine {
namespace MolecularMechanics {
/**
 * @struct AngleType AngleType.h
 * @brief Describes a bond uniquely for given atom types. (useful for maps/unordered_maps)
 *
 *        It is made sure that:
 *           - atoms are ordered correctly.
 *           - two angle types can be compared
 */
struct AngleType {
  AngleType(std::string atom1, std::string atom2, std::string atom3);
  bool operator<(const AngleType& rhs) const;
  bool operator==(const AngleType& rhs) const;

  std::string a1, a2, a3;
};

inline AngleType::AngleType(std::string atom1, std::string atom2, std::string atom3)
  : a1(std::move(atom1)), a2(std::move(atom2)), a3(std::move(atom3)) {
  if (a1 > a3)
    std::swap(a1, a3);
}

inline bool AngleType::operator<(const AngleType& rhs) const {
  if (a1 == rhs.a1) {
    if (a2 == rhs.a2)
      return a3 < rhs.a3;
    else
      return a2 < rhs.a2;
  }
  else
    return a1 < rhs.a1;
}

inline bool AngleType::operator==(const AngleType& rhs) const {
  if ((a1 == rhs.a1) && (a2 == rhs.a2) && (a3 == rhs.a3)) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_ANGLETYPE_H