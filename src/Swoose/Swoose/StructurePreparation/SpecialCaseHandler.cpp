/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SpecialCaseHandler.h"
#include "StructurePreparationData.h"
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace StructurePreparation {
struct StructurePreparationData;

namespace SpecialCaseHandler {

bool isProteinAtom(const StructurePreparationData& data, int index) {
  return std::find(data.vectorOfProteinIndices.begin(), data.vectorOfProteinIndices.end(), index) !=
         data.vectorOfProteinIndices.end();
}

bool isNegative(const StructurePreparationData& data, int index, std::vector<int>& listOfNegatives) {
  if (std::find(listOfNegatives.begin(), listOfNegatives.end(), index) != listOfNegatives.end())
    return false;
  Utils::ElementType element = data.fullStructure.getElement(index);
  auto bondedAtoms = data.listsOfNeighbors[index];
  // Detect carboxylate groups
  if (element == Utils::ElementType::O && bondedAtoms.size() == 1) {
    int neighborIndex = bondedAtoms.front();
    auto neighbor = data.fullStructure.getElement(neighborIndex);
    auto neighborsOfNeighbor = data.listsOfNeighbors[neighborIndex];
    if (neighbor == Utils::ElementType::C && neighborsOfNeighbor.size() == 3) {
      int oxygenCounter = 0;
      for (const auto& neighborOfNeighbor : neighborsOfNeighbor) {
        if (data.fullStructure.getElement(neighborOfNeighbor) == Utils::ElementType::O &&
            data.listsOfNeighbors[neighborOfNeighbor].size() == 1)
          oxygenCounter++;
        if (std::find(listOfNegatives.begin(), listOfNegatives.end(), neighborOfNeighbor) != listOfNegatives.end())
          return false;
      }
      if (oxygenCounter == 2) {
        listOfNegatives.push_back(index);
        return true;
      }
    }
  }
  return false;
}

bool isPositive(const StructurePreparationData& data, int index, std::vector<int>& listOfPositives) {
  if (std::find(listOfPositives.begin(), listOfPositives.end(), index) != listOfPositives.end())
    return false;
  Utils::ElementType element = data.fullStructure.getElement(index);
  auto bondedAtoms = data.listsOfNeighbors[index];
  if ((element == Utils::ElementType::C) && (bondedAtoms.size() == 3)) {
    int nitrogenCounter = 0;
    int nitrogensThatHaveThreeNeighbors = 0;
    for (const auto& neighbor : bondedAtoms) {
      if (data.fullStructure.getElement(neighbor) == Utils::ElementType::N) {
        nitrogenCounter++;
        if (data.listsOfNeighbors[neighbor].size() == 3)
          nitrogensThatHaveThreeNeighbors++;
      }
    }
    if ((nitrogenCounter == 3) && (nitrogensThatHaveThreeNeighbors == 3)) {
      listOfPositives.push_back(index);
      return true;
    }
  }
  // R-NH3+ groups
  else if ((element == Utils::ElementType::N) && (bondedAtoms.size() == 4)) {
    listOfPositives.push_back(index);
    return true;
  }
  return false;
}

bool isCarboxylateC(const StructurePreparationData& data, int index) {
  Utils::ElementType element = data.fullStructure.getElement(index);
  auto neighbors = data.listsOfNeighbors[index];
  if (element == Utils::ElementType::C && neighbors.size() == 3) {
    int oxygenCounter = 0;
    for (const auto& neighbor : neighbors) {
      if (data.fullStructure.getElement(neighbor) == Utils::ElementType::O && data.listsOfNeighbors[neighbor].size() == 1)
        oxygenCounter++;
    }
    if (oxygenCounter == 2)
      return true;
  }
  return false;
}

bool isCTerminus(const StructurePreparationData& data, int index) {
  if (isProteinAtom(data, index) &&
      ((std::find(data.listOfBackboneAlphaCarbons.begin(), data.listOfBackboneAlphaCarbons.end(), index)) !=
       data.listOfBackboneAlphaCarbons.end())) {
    return isCarboxylateC(data, index);
  }
  return false;
}

} // namespace SpecialCaseHandler
} // namespace StructurePreparation
} // namespace Scine
