/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "IndexedStructuralTopologyCreator.h"
#include "IndexedStructuralTopology.h"
#include <Utils/Constants.h>
#include <iostream>

namespace Scine {
namespace MolecularMechanics {

IndexedStructuralTopologyCreator::IndexedStructuralTopologyCreator(const std::vector<std::list<int>>& listsOfNeighbors)
  : nAtoms_(static_cast<int>(listsOfNeighbors.size())), listsOfNeighbors_(listsOfNeighbors) {
}

IndexedStructuralTopology IndexedStructuralTopologyCreator::calculateIndexedStructuralTopology() const {
  IndexedStructuralTopology topology;

  std::set<std::pair<int, int>> excludedNB;
  std::set<std::pair<int, int>> scaledNB;

  for (int a1 = 0; a1 < nAtoms_; ++a1) {
    for (auto a2 : listsOfNeighbors_[a1]) {
      addBond(topology, a1, a2, excludedNB);
      for (auto a3 : listsOfNeighbors_[a2]) {
        if (a3 == a1)
          continue;
        addAngle(topology, a1, a2, a3, excludedNB);
        for (auto a4 : listsOfNeighbors_[a3]) {
          if (a4 == a2)
            continue;
          addDihedral(topology, a1, a2, a3, a4, scaledNB);
        }
        for (auto a4 : listsOfNeighbors_[a2]) {
          if (a4 == a1 || a4 == a3)
            continue;
          if (listsOfNeighbors_[a2].size() == 3)
            addImproperDihedral(topology, a2, a1, a3, a4);
        }
      }
    }
  }

  addNonBondedExclusions(topology, excludedNB, scaledNB);

  return topology;
}

void IndexedStructuralTopologyCreator::addBond(IndexedStructuralTopology& topo, int a1, int a2,
                                               std::set<std::pair<int, int>>& excludedNB) const {
  if (a2 < a1)
    return; // avoids duplication

  topo.addBond(a1, a2);
  excludedNB.emplace(a1, a2);
}

void IndexedStructuralTopologyCreator::addAngle(IndexedStructuralTopology& topo, int a1, int a2, int a3,
                                                std::set<std::pair<int, int>>& excludedNB) const {
  if (a3 < a1)
    return; // avoids duplication

  topo.addAngle(a1, a2, a3);
  excludedNB.emplace(a1, a3);
}

void IndexedStructuralTopologyCreator::addDihedral(IndexedStructuralTopology& topo, int a1, int a2, int a3, int a4,
                                                   std::set<std::pair<int, int>>& scaledNB) const {
  if (a4 < a1)
    return; // avoids duplication

  topo.addDihedral(a1, a2, a3, a4);
  scaledNB.emplace(a1, a4);
}

void IndexedStructuralTopologyCreator::addImproperDihedral(IndexedStructuralTopology& topo, int central, int a2, int a3,
                                                           int a4) const {
  if (a3 < a2 || a4 < a3)
    return; // avoids duplication

  topo.addImproperDihedral(central, a2, a3, a4);
}

void IndexedStructuralTopologyCreator::addNonBondedExclusions(IndexedStructuralTopology& topo,
                                                              const std::set<std::pair<int, int>>& excludedNB,
                                                              std::set<std::pair<int, int>>& scaledNB) const {
  // Verify that non-bonded exclusions have not been counted twice, one time as scaled and one time as fully excluded.
  // In particular, make sure that a the scaled exclusion is erased if it is already a full exclusion.
  // This is needed when cycles are present in the system.

  for (auto p : excludedNB) {
    scaledNB.erase(p);
  }

  // Add to the topology
  for (auto p : excludedNB) {
    topo.addExcludedNonBonded(p.first, p.second);
  }
  for (auto p : scaledNB) {
    topo.addScaledNonBonded(p.first, p.second);
  }
}

void IndexedStructuralTopologyCreator::addHydrogenBondsToIndexedStructuralTopology(IndexedStructuralTopology& topology,
                                                                                   const Utils::AtomCollection& structure) const {
  constexpr double distanceThreshold = 6.0 * Utils::Constants::bohr_per_angstrom; // 6.0 Angstrom
  const auto& elementTypes = structure.getElements();
  std::vector<Utils::ElementType> vectorOfDonorOrAcceptorElements = {Utils::ElementType::N, Utils::ElementType::O,
                                                                     Utils::ElementType::F, Utils::ElementType::Cl};

  const auto& excludedNB = topology.getExcludedNonBondedContainer();
  //  const auto& scaledNB = topology.getScaledNonBondedContainer(); // TODO: see todo below

  // Loop over all bonds
  for (const auto& bond : topology.getBondContainer()) {
    // Loop over possible donor elements
    for (const auto& donorElement : vectorOfDonorOrAcceptorElements) {
      const auto& elementOne = elementTypes[bond.atom1];
      const auto& elementTwo = elementTypes[bond.atom2];
      // Check whether this bond is a bond between H and the donor element
      if ((elementOne == Utils::ElementType::H && elementTwo == donorElement) ||
          (elementTwo == Utils::ElementType::H && elementOne == donorElement)) {
        // Check which atom is the hydrogen
        int hydrogen = bond.atom2;
        int donor = bond.atom1;
        if (elementOne == Utils::ElementType::H) {
          hydrogen = bond.atom1;
          donor = bond.atom2;
        }
        // Loop over all elements
        for (int acceptor = 0; acceptor < elementTypes.size(); ++acceptor) {
          // First, check distance threshold
          if ((structure.getPosition(acceptor) - structure.getPosition(hydrogen)).norm() > distanceThreshold)
            continue;
          if (std::find(vectorOfDonorOrAcceptorElements.begin(), vectorOfDonorOrAcceptorElements.end(),
                        elementTypes[acceptor]) != vectorOfDonorOrAcceptorElements.end()) {
            // Continue if acceptor equals the hydrogen atom
            if (acceptor == hydrogen)
              continue;
            // Continue if acceptor equals the donor atom
            if (acceptor == donor)
              continue;
            // Continue if acceptor and donor are in the excluded non-bonded container (1,2- or 1,3-connection)
            if (std::find(excludedNB.begin(), excludedNB.end(), StructuralExcludedNonBonded(donor, acceptor)) !=
                excludedNB.end())
              continue;
            // TODO: Continue if acceptor and donor are in the scaled non-bonded container (1,4-connection) ???
            /*
            if (std::find(scaledNB.begin(), scaledNB.end(), IndexedStructuralScaledNonBonded(donor, acceptor)) !=
                scaledNB.end())
              continue;
            */
            addHydrogenBond(topology, donor, hydrogen, acceptor);
          }
        }
      }
    }
  }
}

void IndexedStructuralTopologyCreator::addHydrogenBond(IndexedStructuralTopology& topo, int a1, int a2, int a3) const {
  topo.addHydrogenBond(a1, a2, a3);
}

} // namespace MolecularMechanics
} // namespace Scine