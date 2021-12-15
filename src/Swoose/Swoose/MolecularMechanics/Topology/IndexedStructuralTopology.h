/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_INDEXEDSTRUCTURALTOPOLOGY_H
#define MOLECULARMECHANICS_INDEXEDSTRUCTURALTOPOLOGY_H

#include "IndexedStructuralElements.h"
#include <vector>

namespace Scine {
namespace MolecularMechanics {

/**
 * @class IndexedStructuralTopology IndexedStructuralTopology.h
 * @brief Class containing the structural information about the connectivity of a system.
 *
 *        It holds lists of bonds, angles, dihedrals, etc. for MM calculations.
 *        It does not contain information about atom types
 */

class IndexedStructuralTopology {
 public:
  /** @brief Adds a bond. Performs no check whether it already exists. */
  void addBond(int a1, int a2);
  /** @brief Adds a angle. Performs no check whether it already exists. */
  void addAngle(int a1, int a2, int a3);
  /** @brief Adds a dihedral. Performs no check whether it already exists. */
  void addDihedral(int a1, int a2, int a3, int a4);
  /** @brief Adds an improper dihedral. Performs no check whether it already exists. */
  void addImproperDihedral(int central, int a2, int a3, int a4);
  /** @brief Adds an excluded non-bonded interaction. Performs no check whether it already exists. */
  void addExcludedNonBonded(int a1, int a2);
  /** @brief Adds a scaled excluded non-bonded interaction. Performs no check whether it already exists. */
  void addScaledNonBonded(int a1, int a2);
  /** @brief Adds a scaled excluded non-bonded interaction. Performs no check whether it already exists. */
  void addHydrogenBond(int a1, int a2, int a3);

  /** @name Members returning the arrays of structural elements
   * @{
   */
  const std::vector<IndexedStructuralBond>& getBondContainer() const;
  std::vector<IndexedStructuralBond>& getBondContainer();
  const std::vector<IndexedStructuralAngle>& getAngleContainer() const;
  std::vector<IndexedStructuralAngle>& getAngleContainer();
  const std::vector<IndexedStructuralDihedral>& getDihedralContainer() const;
  std::vector<IndexedStructuralDihedral>& getDihedralContainer();
  const std::vector<IndexedStructuralImproperDihedral>& getImproperDihedralContainer() const;
  std::vector<IndexedStructuralImproperDihedral>& getImproperDihedralContainer();
  const std::vector<StructuralExcludedNonBonded>& getExcludedNonBondedContainer() const;
  const std::vector<IndexedStructuralScaledNonBonded>& getScaledNonBondedContainer() const;
  const std::vector<IndexedStructuralHydrogenBond>& getHydrogenBondContainer() const;
  /*! @} */

 private:
  std::vector<IndexedStructuralBond> bonds_;
  std::vector<IndexedStructuralAngle> angles_;
  std::vector<IndexedStructuralDihedral> dihedrals_;
  std::vector<IndexedStructuralImproperDihedral> improperDihedrals_;
  std::vector<StructuralExcludedNonBonded> excludedNB_;
  std::vector<IndexedStructuralScaledNonBonded> scaledNB_;
  std::vector<IndexedStructuralHydrogenBond> hydrogenBonds_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_INDEXEDSTRUCTURALTOPOLOGY_H