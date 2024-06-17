/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_TOPOLOGYUTILS_H
#define SWOOSEUTILITIES_TOPOLOGYUTILS_H

#include <Utils/Typenames.h>
#include <list>
#include <vector>

namespace Scine {

namespace Utils {
class BondOrderCollection;
class Atom;
class AtomCollection;
} // namespace Utils

namespace SwooseUtilities {

class TopologyUtils {
 public:
  /**
   * @brief This function generates an array of nAtoms lists of neighbors from a lower-triangular bond-order matrix.
   */
  static std::vector<std::list<int>>
  generateListsOfNeighborsFromBondOrderMatrix(int nAtoms, const Utils::BondOrderCollection& bondOrderMatrix,
                                              double minimalBondOrderToConsider);
  /**
   * @brief This function generates a bond order matrix from a given lists of neighbors.
   */
  static Utils::BondOrderCollection generateBondOrderMatrixFromListsOfNeighbors(const std::vector<std::list<int>>& listsOfNeighbors);
  /**
   * @brief Divides the atoms into two groups according to which side of the bond given by atoms a1 and a2 are.
   *
   *        In the returned array, index i is 1 if atom i is on the side of a1, 2 if it is on the side of a2, and 0 if
   *        if elsewhere. ok is true if the division was successful (i.e. no cycles).
   */
  static std::vector<int> divideStructureAtBond(int a1, int a2, const Utils::BondOrderCollection& bondOrderMatrix, bool& ok);
  /**
   * @brief Counts and returns the number of neighbors for atom with index 'atomIndex' that has the given element type.
   */
  static int countNeighborsOfElementType(int atomType, const std::vector<std::list<int>>& listsOfNeighbors,
                                         Utils::ElementType elementType,
                                         const Utils::ElementTypeCollection& elementTypeCollection);
  /**
   * @brief
   *
   * @param nAtoms
   * @param listsOfNeighbors
   * @param nNeighbors
   */
  static void calculateNumberOfNeighbors(int nAtoms, std::vector<std::list<int>> listsOfNeighbors, std::vector<int>& nNeighbors);
};

} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_TOPOLOGYUTILS_H
