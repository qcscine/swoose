/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FragmentationHelper.h"
#include <Core/Log.h>
#include <Molassembler/Graph.h>
#include <Molassembler/Interpret.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Geometry/GeometryUtilities.h>

namespace Scine {
namespace SwooseUtilities {
namespace FragmentationHelper {

std::mt19937 randomEngine(42);

void addAtomsUpToReasonableCut(Utils::AtomCollection& atomsToAdd, std::vector<int>& alreadyAddedAtoms,
                               std::deque<bool>& isSaturatingAtom, int atomInside, int atomOutside,
                               const Utils::AtomCollection& fullStructure, const std::vector<std::list<int>>& listsOfNeighbors,
                               double probabilityToDivide, std::shared_ptr<std::mt19937> randomEngine) {
  // If both atoms are already in the atomsToAdd structure, just return from this function.
  // This prevents problems with cycles.
  if (std::find(alreadyAddedAtoms.begin(), alreadyAddedAtoms.end(), atomInside) != alreadyAddedAtoms.end()) {
    if (std::find(alreadyAddedAtoms.begin(), alreadyAddedAtoms.end(), atomOutside) != alreadyAddedAtoms.end()) {
      return;
    }
  }
  else {
    alreadyAddedAtoms.push_back(atomInside); // add inside atom to already added atoms vector
  }

  // First, check whether the system is divisible at this bond
  if (!isDivisibleAtBond(atomInside, atomOutside, fullStructure, listsOfNeighbors, probabilityToDivide, randomEngine)) {
    // Add the outside atom to the temporary atomsToAdd structure
    atomsToAdd.push_back(fullStructure[atomOutside]);
    isSaturatingAtom.push_back(false);        // The added atom is not an atom placed for saturation
    alreadyAddedAtoms.push_back(atomOutside); // add outside atom to already added atoms vector

    // Loop over all neighbors of that outside atom
    for (const auto& neighbor : listsOfNeighbors[atomOutside]) {
      if (neighbor == atomInside)
        continue;
      // Call the function again making the outside atom the new inside atom and the neighbor the outside atom
      addAtomsUpToReasonableCut(atomsToAdd, alreadyAddedAtoms, isSaturatingAtom, atomOutside, neighbor, fullStructure,
                                listsOfNeighbors, probabilityToDivide, randomEngine);
    }
  }
  else {
    // If it is, cut at that bond and saturate with hydrogen
    Utils::ElementType addedElement = Utils::ElementType::H;
    Eigen::RowVector3d bondVector = fullStructure.getPosition(atomOutside) - fullStructure.getPosition(atomInside);

    // Covalent radii
    Utils::ElementType elementAtomInside = fullStructure[atomInside].getElementType();
    auto covalentRadius1 = Utils::ElementInfo::covalentRadius(elementAtomInside);
    auto covalentRadius2 = Utils::ElementInfo::covalentRadius(addedElement);

    // add replacement atom along the bond vector of the bond "insideAtom - outsideAtom"
    Utils::Atom addedAtom;
    Eigen::RowVector3d positionOfAtomInside = fullStructure.getPosition(atomInside);

    Eigen::RowVector3d scaledBondVector = ((covalentRadius1 + covalentRadius2) / bondVector.norm()) * bondVector;
    Eigen::RowVector3d addedAtomPosition = positionOfAtomInside + scaledBondVector;

    addedAtom.setElementType(addedElement);
    addedAtom.setPosition(addedAtomPosition);
    atomsToAdd.push_back(addedAtom);
    isSaturatingAtom.push_back(true); // The added atom is an atom placed for saturation
  }
}

void addMoreAtomsToSubsystem(const Utils::AtomCollection& atomsToAdd, Utils::AtomCollection& subsystem,
                             const std::deque<bool>& isSaturatingAtom, bool addSaturatingAtoms) {
  int atomCounter = 0;
  for (const auto& atom : atomsToAdd) {
    // Check whether atom is an atom placed for valence saturation.
    bool thisAtomIsSaturatingAtom = isSaturatingAtom[atomCounter];

    // Add saturating atoms
    if (thisAtomIsSaturatingAtom && addSaturatingAtoms) {
      // Check that the atom that is about to be added is not close to an existing atom
      if (!atomIsCloseToExistingAtom(atom, subsystem)) {
        subsystem.push_back(atom);
      }
    }

    // Add non-saturating atoms
    else if ((!thisAtomIsSaturatingAtom) && (!addSaturatingAtoms)) {
      // Check that the atom that is about to be added is not close to an existing atom
      if (!atomIsCloseToExistingAtom(atom, subsystem)) {
        subsystem.push_back(atom);
      }
    }
    atomCounter++;
  }
}

void checkSizeOfSubsystem(int size, double& additionToRadius, bool& unsuccessfulFragmentation, int fragmentIndex,
                          int fullStructureSize, int maximumSubsystemSize, Core::Log& log) {
  if (!unsuccessfulFragmentation) {
    if (size < minimumSubsystemSize_ && size != fullStructureSize) {
      additionToRadius += 1.0; // add 1 bohr to radius if subsystem is very small
      unsuccessfulFragmentation = true;
    }
    else if (size > maximumSubsystemSize) {
      log.warning << "Size of fragment centered around atom " << fragmentIndex
                  << " is perhaps too large. Number of atoms: " << size << Core::Log::endl;
    }
  }
}

bool atomIsCloseToExistingAtom(const Utils::Atom& atom, const Utils::AtomCollection& existingAtoms) {
  constexpr double threshold = 0.7 * Utils::Constants::bohr_per_angstrom; // TODO: Find optimal value
  bool result = false;
  for (const auto& a : existingAtoms) {
    double distance = (atom.getPosition() - a.getPosition()).norm();
    if (distance < threshold) {
      result = true;
      break;
    }
  }
  return result;
}

// TODO: Improve this function
bool isDivisibleAtBond(int atomInside, int atomOutside, const Utils::AtomCollection& fullStructure,
                       const std::vector<std::list<int>>& listsOfNeighbors, double probabilityToDivide,
                       std::shared_ptr<std::mt19937> randomEngine) {
  assert(randomEngine != nullptr && "Random engine was not provided.");
  bool ret = false;
  Utils::ElementType elementInside = fullStructure[atomInside].getElementType();
  Utils::ElementType elementOutside = fullStructure[atomOutside].getElementType();

  if (elementInside == Utils::ElementType::C && elementOutside == Utils::ElementType::C)
    ret = listsOfNeighbors[atomInside].size() == 4;
  else if (elementInside == Utils::ElementType::C && elementOutside == Utils::ElementType::N)
    ret = listsOfNeighbors[atomInside].size() == 4;

  if (ret && probabilityToDivide < 1.0) {
    std::bernoulli_distribution distribution(probabilityToDivide);
    return distribution(*randomEngine);
  }

  return ret;
}

void updateInformationForIndexMapping(const Utils::AtomCollection& subsystem,
                                      const Utils::AtomCollection& fullStructure, std::vector<int>& atomIndexMapping) {
  std::vector<int> indicesOfAtomsInFragment;

  int atomIndex = 0;
  for (const auto& atom : subsystem) {
    try {
      // Determine the index of the fragment's atom in the full system
      auto indexOfAtomInFullSystem = Utils::Geometry::Distances::getIndexOfAtomInStructure(fullStructure, atom);
      indicesOfAtomsInFragment.push_back(indexOfAtomInFullSystem);
    }
    // The above function throws a runtime_error when the atom is not found in the given full system's structure
    catch (const std::runtime_error& e) {
      // -1 is the index of all atoms that are not in the full structure, but in the fragment
      indicesOfAtomsInFragment.push_back(indexForSaturatingAtoms);
    }
    atomIndex++;
  }
  atomIndexMapping = indicesOfAtomsInFragment;
}

std::vector<int> calculateSubgraphSizes(const Utils::AtomCollection& fullStructure, const Utils::BondOrderCollection& bondOrders) {
  std::vector<int> subgraphSizes;
  auto numberOfAtoms = fullStructure.size();
  auto bondInterpretation = Molassembler::Interpret::graphs(fullStructure, bondOrders);
  subgraphSizes.reserve(numberOfAtoms);
  for (int i = 0; i < numberOfAtoms; ++i) {
    auto subgraphSize = bondInterpretation.graphs[bondInterpretation.componentMap.apply(i).component].V();
    subgraphSizes.push_back(subgraphSize);
  }
  return subgraphSizes;
}

} // namespace FragmentationHelper
} // namespace SwooseUtilities
} // namespace Scine
