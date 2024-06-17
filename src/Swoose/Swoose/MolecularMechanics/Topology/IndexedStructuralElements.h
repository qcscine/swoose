/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_INDEXEDSTRUCTURALELEMENTS_H
#define MOLECULARMECHANICS_INDEXEDSTRUCTURALELEMENTS_H

namespace Scine {
namespace MolecularMechanics {

/**
 * @brief Holds indexes for a bond.
 */
struct IndexedStructuralBond {
  IndexedStructuralBond(int a1, int a2) : atom1(a1), atom2(a2) {
  }

  int atom1;
  int atom2;
};

/**
 * @brief Holds indexes for an angle.
 *        The connectivity is atom1 -- atom2 -- atom3
 */
struct IndexedStructuralAngle {
  IndexedStructuralAngle(int a1, int a2, int a3) : atom1(a1), atom2(a2), atom3(a3) {
  }

  int atom1;
  int atom2;
  int atom3;
};

/**
 * @brief Holds indexes for a dihedral.
 *        The connectivity is atom1 -- atom2 -- atom3 -- atom4
 */
struct IndexedStructuralDihedral {
  IndexedStructuralDihedral(int a1, int a2, int a3, int a4) : atom1(a1), atom2(a2), atom3(a3), atom4(a4) {
  }

  int atom1;
  int atom2;
  int atom3;
  int atom4;

  bool operator==(const IndexedStructuralDihedral& other) const {
    return ((atom1 == other.atom1) && (atom2 == other.atom2) && (atom3 == other.atom3) && (atom4 == other.atom4));
  }
};

/**
 * @brief Holds indexes for an improper dihedral.
 */
struct IndexedStructuralImproperDihedral {
  IndexedStructuralImproperDihedral(int central, int a2, int a3, int a4)
    : centralAtom(central), atom2(a2), atom3(a3), atom4(a4) {
  }

  int centralAtom;
  int atom2;
  int atom3;
  int atom4;
};

/**
 * @brief Holds indexes for an excluded non-bonded interaction. (1-2 and 1-3 neighbors)
 */
struct StructuralExcludedNonBonded {
  StructuralExcludedNonBonded(int a1, int a2) : atom1(a1), atom2(a2) {
  }

  int atom1;
  int atom2;

  bool operator==(const StructuralExcludedNonBonded& other) const {
    return (((atom1 == other.atom1) && (atom2 == other.atom2)) || ((atom1 == other.atom2) && (atom2 == other.atom1)));
  }
};

/**
 * @brief Holds indexes for a potentially scaled non-bonded interaction. (1-4 neighbors)
 */
struct IndexedStructuralScaledNonBonded {
  IndexedStructuralScaledNonBonded(int a1, int a2) : atom1(a1), atom2(a2) {
  }

  int atom1;
  int atom2;

  bool operator==(const IndexedStructuralScaledNonBonded& other) const {
    return (((atom1 == other.atom1) && (atom2 == other.atom2)) || ((atom1 == other.atom2) && (atom2 == other.atom1)));
  }
};

/**
 * @brief Holds indexes for a hydrogen-bond interaction.
 */
struct IndexedStructuralHydrogenBond {
  IndexedStructuralHydrogenBond(int a1, int a2, int a3) : atom1(a1), atom2(a2), atom3(a3) {
  }

  int atom1;
  int atom2;
  int atom3;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_INDEXEDSTRUCTURALELEMENTS_H
