/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_INDEXEDSTRUCTURALTOPOLOGYCREATOR_H
#define MOLECULARMECHANICS_INDEXEDSTRUCTURALTOPOLOGYCREATOR_H

#include <Utils/Geometry/AtomCollection.h>
#include <list>
#include <set>
#include <vector>

namespace Scine {

namespace MolecularMechanics {
class IndexedStructuralTopology;

/**
 * @class IndexedStructuralTopologyCreator IndexedStructuralTopologyCreator.h
 * @brief Class that creates an IndexedStructuralTopology based on the connectivity of a molecular system.
 *
 *        It finds bonds, angles, dihedrals, etc. in a molecular system. Terms are added with no duplication.
 *        It does not determine atom types.
 */
class IndexedStructuralTopologyCreator {
 public:
  /** @brief Constructor taking as parameter a list of neighbors for each atom. */
  explicit IndexedStructuralTopologyCreator(const std::vector<std::list<int>>& listsOfNeighbors);

  /** @brief Function generating an IndexedStructuralTopology and returning it. */
  IndexedStructuralTopology calculateIndexedStructuralTopology() const;
  /** @brief Function adding hydrogen bonds to the IndexedStructuralTopology that is given as an argument */
  void addHydrogenBondsToIndexedStructuralTopology(IndexedStructuralTopology& topology,
                                                   const Utils::AtomCollection& structure) const;

 private:
  void addBond(IndexedStructuralTopology& topo, int a1, int a2, std::set<std::pair<int, int>>& excludedNB) const;
  void addAngle(IndexedStructuralTopology& topo, int a1, int a2, int a3, std::set<std::pair<int, int>>& excludedNB) const;
  void addDihedral(IndexedStructuralTopology& topo, int a1, int a2, int a3, int a4, std::set<std::pair<int, int>>& scaledNB) const;
  void addImproperDihedral(IndexedStructuralTopology& topo, int central, int a2, int a3, int a4) const;
  void addNonBondedExclusions(IndexedStructuralTopology& topo, const std::set<std::pair<int, int>>& excludedNB,
                              std::set<std::pair<int, int>>& scaledNB) const;
  void addHydrogenBond(IndexedStructuralTopology& topo, int a1, int a2, int a3) const;

  int nAtoms_;
  std::vector<std::list<int>> listsOfNeighbors_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_INDEXEDSTRUCTURALTOPOLOGYCREATOR_H
