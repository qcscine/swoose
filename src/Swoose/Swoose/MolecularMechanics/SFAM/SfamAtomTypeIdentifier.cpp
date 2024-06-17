/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SfamAtomTypeIdentifier.h"
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace MolecularMechanics {

SfamAtomTypeIdentifier::SfamAtomTypeIdentifier(int nAtoms, Utils::ElementTypeCollection elementTypes,
                                               const std::vector<std::list<int>>& listsOfNeighbors)
  : nAtoms_(nAtoms), elementTypes_(std::move(elementTypes)), listsOfNeighbors_(std::move(listsOfNeighbors)) {
  atomTypes_.resize(nAtoms);
  assert((nAtoms_ == elementTypes_.size() && nAtoms_ == static_cast<int>(listsOfNeighbors_.size())) &&
         "Arrays have wrong size in atom type identifier (not equal to the number of atoms).");
}

void SfamAtomTypeIdentifier::setAtomType(int index, std::string type) {
  atomTypes_[index] = std::move(type);
}

AtomTypesHolder SfamAtomTypeIdentifier::getAtomTypes(SfamAtomTypeLevel atl) {
  if (atl == SfamAtomTypeLevel::Low) {
    for (int i = 0; i < nAtoms_; ++i) {
      auto atomSymbol = Utils::ElementInfo::symbol(elementTypes_[i]);
      auto numberOfNeighbors = listsOfNeighbors_[i].size();
      setAtomType(i, atomSymbol + std::to_string(numberOfNeighbors));
    }
  }

  if (atl == SfamAtomTypeLevel::High) {
    for (int i = 0; i < nAtoms_; ++i) {
      auto atomSymbol = Utils::ElementInfo::symbol(elementTypes_[i]);
      auto numberOfNeighbors = listsOfNeighbors_[i].size();
      if (atomSymbol == "H") {
        setAtomType(i, "H");
      }
      else {
        std::vector<Utils::ElementType> neighborAtoms;
        for (const int& neighbor : listsOfNeighbors_[i]) {
          auto elementTypeOfNeighbor = elementTypes_[neighbor];
          if (std::find(neighborAtoms.begin(), neighborAtoms.end(), elementTypeOfNeighbor) == neighborAtoms.end()) {
            neighborAtoms.push_back(elementTypeOfNeighbor);
          }
        }
        std::sort(neighborAtoms.begin(), neighborAtoms.end());
        std::string atomType = atomSymbol + std::to_string(numberOfNeighbors) + "_";
        for (const auto& neighbor : neighborAtoms) {
          int nNeighborsOfElementType =
              SwooseUtilities::TopologyUtils::countNeighborsOfElementType(i, listsOfNeighbors_, neighbor, elementTypes_);
          atomType += Utils::ElementInfo::symbol(neighbor) + std::to_string(nNeighborsOfElementType);
        }
        setAtomType(i, atomType);
      }
    }

    for (int j = 0; j < int(atomTypes_.size()); ++j) {
      if (atomTypes_[j] == "H") {
        std::string atomType = "H_";
        for (const auto& neighbor : listsOfNeighbors_[j]) {
          atomType += atomTypes_[neighbor];
        }
        setAtomType(j, atomType);
      }
    }
  }

  if (atl == SfamAtomTypeLevel::Elements) {
    for (int i = 0; i < nAtoms_; ++i) {
      auto atomSymbol = Utils::ElementInfo::symbol(elementTypes_[i]);
      setAtomType(i, atomSymbol);
    }
  }

  if (atl == SfamAtomTypeLevel::Unique) {
    for (int i = 0; i < nAtoms_; ++i) {
      auto atomSymbol = Utils::ElementInfo::symbol(elementTypes_[i]);
      setAtomType(i, atomSymbol + std::to_string(i));
    }
  }

  return AtomTypesHolder(atomTypes_);
}

SfamAtomTypeLevel SfamAtomTypeIdentifier::generateSfamAtomTypeLevelFromString(std::string sfamAtomTypeLevelString) {
  if (sfamAtomTypeLevelString == "elements")
    return SfamAtomTypeLevel::Elements;
  else if (sfamAtomTypeLevelString == "low")
    return SfamAtomTypeLevel::Low;
  else if (sfamAtomTypeLevelString == "high")
    return SfamAtomTypeLevel::High;
  else if (sfamAtomTypeLevelString == "unique")
    return SfamAtomTypeLevel::Unique;
  else
    throw std::runtime_error("You did not specify a valid atom type level.");
}

} // namespace MolecularMechanics
} // namespace Scine
