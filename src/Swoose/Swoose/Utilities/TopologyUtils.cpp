/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TopologyUtils.h"
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/Geometry/AtomCollection.h>
#include <queue>

namespace Scine {
namespace SwooseUtilities {

std::vector<std::list<int>>
TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(int nAtoms, const Utils::BondOrderCollection& bondOrderMatrix,
                                                           double minimalBondOrderToConsider) {
  assert(nAtoms >= 0);

  std::vector<std::list<int>> listsOfNeighbors(nAtoms);

  for (int i = 1; i < nAtoms; ++i) {
    for (int j = 0; j < i; ++j) {
      double bondOrder = bondOrderMatrix.getOrder(i, j);
      if (bondOrder > minimalBondOrderToConsider) {
        listsOfNeighbors[i].push_back(j);
        listsOfNeighbors[j].push_back(i);
      }
    }
  }

  return listsOfNeighbors;
}

Utils::BondOrderCollection
TopologyUtils::generateBondOrderMatrixFromListsOfNeighbors(const std::vector<std::list<int>>& listsOfNeighbors) {
  Utils::BondOrderCollection bondOrderMatrix(listsOfNeighbors.size());
  for (int i = 0; i < listsOfNeighbors.size(); ++i)
    for (const auto& j : listsOfNeighbors[i])
      bondOrderMatrix.setOrder(i, j, 1.0);

  return bondOrderMatrix;
}

std::vector<int> TopologyUtils::divideStructureAtBond(int a1, int a2, const Utils::BondOrderCollection& bondOrderMatrix,
                                                      bool& ok) {
  auto nAtoms = bondOrderMatrix.getSystemSize();
  std::vector<int> classes(nAtoms, 0);
  classes[a1] = 1;
  classes[a2] = 2;

  auto listsOfNeighbors = generateListsOfNeighborsFromBondOrderMatrix(nAtoms, bondOrderMatrix, 0.5);

  //                         index, parent, group
  using Element = std::tuple<int, int, int>;
  std::queue<Element> q;

  // Add first elements: neighbors of a1 and a2
  for (auto n : listsOfNeighbors[a1])
    if (n != a2)
      q.emplace(n, a1, 1);
  for (auto n : listsOfNeighbors[a2])
    if (n != a1)
      q.emplace(n, a2, 2);

  while (!q.empty()) {
    auto e = q.front();
    q.pop();

    auto index = std::get<0>(e);
    auto parent = std::get<1>(e);
    auto group = std::get<2>(e);

    // Check if it was already marked
    if (classes[index] == group)
      continue;

    // Check that it is not a cycle
    if (classes[index] != 0) {
      ok = false;
      return classes;
    }

    classes[index] = group;
    for (auto n : listsOfNeighbors[index])
      if (n != parent)
        q.emplace(n, index, group);
  }

  ok = true;
  return classes;
}

int TopologyUtils::countNeighborsOfElementType(int atomType, const std::vector<std::list<int>>& listsOfNeighbors,
                                               Utils::ElementType elementType,
                                               const Utils::ElementTypeCollection& elementTypeCollection) {
  if (atomType >= listsOfNeighbors.size() || listsOfNeighbors.size() != elementTypeCollection.size())
    throw std::runtime_error("Incompatible function arguments for counting neighbors of given element type.");
  int nNeighborsOfType = 0;
  for (int n : listsOfNeighbors[atomType]) {
    if (elementTypeCollection.at(n) == elementType)
      nNeighborsOfType++;
  }
  return nNeighborsOfType;
}

} // namespace SwooseUtilities
} // namespace Scine