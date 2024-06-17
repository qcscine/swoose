/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaffAtomTypeIdentifier.h"
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Geometry/ElementInfo.h>
#include <fstream>
#include <stack>
#include <tuple>
#include <utility>

namespace Scine {
namespace MolecularMechanics {

GaffAtomTypeIdentifier::GaffAtomTypeIdentifier(int nAtoms, Utils::ElementTypeCollection elementTypes,
                                               const std::vector<std::list<int>>& listsOfNeighbors, std::string atomTypesFile)
  : nAtoms_(nAtoms),
    elementTypes_(std::move(elementTypes)),
    listsOfNeighbors_(std::move(listsOfNeighbors)),
    nNeighbors_(nAtoms),
    numberOfCyclesForAtom_(nAtoms, 0),
    gaffAtomTypes_(nAtoms),
    isConjugated_(nAtoms, false),
    atomTypesFile_(atomTypesFile) {
  assert((nAtoms_ == elementTypes_.size() && nAtoms_ == static_cast<int>(listsOfNeighbors_.size())) &&
         "Arrays have wrong size in atom type identifier (not equal to the number of atoms).");
}

void GaffAtomTypeIdentifier::setAtomType(int index, std::string type, bool throwIfAlreadySet) {
  if (throwIfAlreadySet && atomTypeSet(index))
    throw std::runtime_error("Atom type of MM atom was already set.");
  gaffAtomTypes_.at(index) = std::move(type);
}

AtomTypesHolder GaffAtomTypeIdentifier::getAtomTypes() {
  // Read from file if such a file is available
  if (!atomTypesFile_.empty()) {
    readGaffAtomTypesFromFile();
    return AtomTypesHolder(gaffAtomTypes_);
  }

  // Automated algorithm (still beta version, not working perfectly)
  setArraysForElementTypes();
  SwooseUtilities::TopologyUtils::calculateNumberOfNeighbors(nAtoms_, listsOfNeighbors_, nNeighbors_);

  handleHalogens();
  handleOxygenFunctionalGroups();
  handleSulfurFunctionalGroups();
  handleSpecialNitrogens();
  handlePhosphorus();
  handleCarbons();

  lookForCycles();
  handleCycles();

  checkConjugation();
  handleRemainingN();
  handleHydrogens();
  return AtomTypesHolder(gaffAtomTypes_);
}

bool GaffAtomTypeIdentifier::atomTypeSet(int atomIndex) const {
  return !gaffAtomTypes_.at(atomIndex).empty();
}

void GaffAtomTypeIdentifier::setArraysForElementTypes() {
  for (int a = 0; a < nAtoms_; ++a) {
    Utils::ElementType element = elementTypes_[a];
    if (element == Utils::ElementType::H)
      atomsH_.push_back(a);
    else if (element == Utils::ElementType::C)
      atomsC_.push_back(a);
    else if (element == Utils::ElementType::N)
      atomsN_.push_back(a);
    else if (element == Utils::ElementType::O)
      atomsO_.push_back(a);
    else if (element == Utils::ElementType::P)
      atomsP_.push_back(a);
    else if (element == Utils::ElementType::S)
      atomsS_.push_back(a);
    else if (element == Utils::ElementType::F)
      atomsF_.push_back(a);
    else if (element == Utils::ElementType::Cl)
      atomsCl_.push_back(a);
    else if (element == Utils::ElementType::Br)
      atomsBr_.push_back(a);
    else if (element == Utils::ElementType::I)
      atomsI_.push_back(a);
    else
      throw std::runtime_error("Element type is not supported by GAFF.");
  }
}

void GaffAtomTypeIdentifier::verifyNeighborNumber(int atomIndex, int imposedNumberOfNeighbors, bool upper_limit) {
  bool throwError = (upper_limit && nNeighbors_[atomIndex] > imposedNumberOfNeighbors) ||
                    (!upper_limit && nNeighbors_[atomIndex] != imposedNumberOfNeighbors);
  if (throwError) {
    throw std::runtime_error("GAFF atom type detection failed! MM atom does not have a suitable number of neighbors."
                             " Please provide the atom types manually! Atom index: " +
                             std::to_string(atomIndex));
  }
}

void GaffAtomTypeIdentifier::handleHalogens() {
  for (auto a : atomsF_) {
    setAtomType(a, "f");
    verifyNeighborNumber(a, 1, true);
  }
  for (auto a : atomsCl_) {
    setAtomType(a, "cl");
    verifyNeighborNumber(a, 1, true);
  }
  for (auto a : atomsBr_) {
    setAtomType(a, "br");
    verifyNeighborNumber(a, 1, true);
  }
  for (auto a : atomsI_) {
    setAtomType(a, "i");
    verifyNeighborNumber(a, 1, true);
  }
}

void GaffAtomTypeIdentifier::handleOxygenFunctionalGroups() {
  for (auto oAtom : atomsO_) {
    if (atomTypeSet(oAtom))
      continue;
    int neighbors = nNeighbors_[oAtom];

    // Check water, hydronium, hydroxyl
    int nHydrogens = SwooseUtilities::TopologyUtils::countNeighborsOfElementType(oAtom, listsOfNeighbors_,
                                                                                 Utils::ElementType::H, elementTypes_);
    if (neighbors == nHydrogens) {
      setAtomType(oAtom, "ow");
      for (auto h : listsOfNeighbors_[oAtom]) {
        if (elementTypes_[h] == Utils::ElementType::H)
          setAtomType(h, "hw");
      }
      continue;
    }

    if (neighbors > 2 || neighbors == 0)
      throw std::runtime_error("MM oxygen atom does not have suitable number of neighbors.");

    if (neighbors == 2) {
      if (nHydrogens == 1) { // i.e. hydroxy
        setAtomType(oAtom, "oh");
        for (auto h : listsOfNeighbors_[oAtom]) {
          if (elementTypes_[h] == Utils::ElementType::H)
            setAtomType(h, "ho");
        }
      }
      else { // i.e. ester, ether, etc
        setAtomType(oAtom, "os");
      }
    }
    else { // i.e. one neighbor
      int nIndex = listsOfNeighbors_[oAtom].front();
      auto nType = elementTypes_[nIndex];

      // check if it's deprotonated hydroxy:
      if (nType == Utils::ElementType::C && nNeighbors_[nIndex] == 4) {
        setAtomType(oAtom, "oh");
      }
      else {
        setAtomType(oAtom, "o");

        // carbonyl
        if (nType == Utils::ElementType::C)
          setAtomType(nIndex, "c", false);
        // nitro
        else if (nType == Utils::ElementType::N)
          setAtomType(nIndex, "no", false);
      }
    }
  }
}

void GaffAtomTypeIdentifier::handleSulfurFunctionalGroups() {
  for (auto sAtom : atomsS_) {
    if (atomTypeSet(sAtom))
      continue;
    int neighbors = nNeighbors_[sAtom];

    if (neighbors == 2) {
      int nHydrogens = SwooseUtilities::TopologyUtils::countNeighborsOfElementType(sAtom, listsOfNeighbors_,
                                                                                   Utils::ElementType::H, elementTypes_);
      if (nHydrogens == 1) {
        setAtomType(sAtom, "sh");
        for (auto h : listsOfNeighbors_[sAtom]) {
          if (elementTypes_[h] == Utils::ElementType::H)
            setAtomType(h, "hs");
        }
      }
      else {
        setAtomType(sAtom, "ss");
      }
    }
    else if (neighbors == 1) {
      setAtomType(sAtom, "s");

      int nIndex = listsOfNeighbors_[sAtom].front();
      auto nType = elementTypes_[nIndex];

      if (nType == Utils::ElementType::C)
        setAtomType(nIndex, "c");
    }
    // TODO: hypervalent, check better!!!
    else if (neighbors == 3)
      setAtomType(sAtom, "s4");
    else if (neighbors == 4)
      setAtomType(sAtom, "s6");
  }
}

void GaffAtomTypeIdentifier::handleSpecialNitrogens() {
  for (auto a : atomsN_) {
    if (atomTypeSet(a))
      continue;

    if (nNeighbors_[a] == 4) {
      setAtomType(a, "n4");
      continue;
    }

    if (nNeighbors_[a] == 2) {
      setAtomType(a, "n2");
      continue;
    }

    bool hasCarbonylLikeNeighbor = false;
    for (auto n : listsOfNeighbors_[a]) {
      if (gaffAtomTypes_[n] == "c")
        hasCarbonylLikeNeighbor = true;
      else if (elementTypes_[n] == Utils::ElementType::S) {
        for (auto n2 : listsOfNeighbors_[n]) {
          if (gaffAtomTypes_[n2] == "o")
            hasCarbonylLikeNeighbor = true;
        }
      }
    }

    if (hasCarbonylLikeNeighbor)
      setAtomType(a, "n");
  }
}

void GaffAtomTypeIdentifier::handlePhosphorus() {
  for (auto a : atomsP_) {
    if (nNeighbors_[a] == 2) {
      setAtomType(a, "p2");
      continue;
    }

    bool doubleBondPossible = false;
    for (auto n : listsOfNeighbors_[a]) {
      auto nType = gaffAtomTypes_[n];
      if (nType == "n2" || nType == "o" || nType == "s2")
        doubleBondPossible = true;
      else if (elementTypes_[n] == Utils::ElementType::C && nNeighbors_[n] < 3 && nType != "c")
        doubleBondPossible = true;
      // TODO: else if P or S ?
    }

    if (!doubleBondPossible && nNeighbors_[a] == 3) {
      setAtomType(a, "p3");
      continue;
    }
    else if (nNeighbors_[a] == 3)
      setAtomType(a, "p4");
    else if (nNeighbors_[a] == 4)
      setAtomType(a, "p5");
  }
}

void GaffAtomTypeIdentifier::handleCarbons() {
  for (auto a : atomsC_) {
    if (atomTypeSet(a))
      continue;

    int neighbors = nNeighbors_[a];

    if (neighbors == 2)
      setAtomType(a, "c1");
    else if (neighbors == 3)
      setAtomType(a, "c2");
    else if (neighbors == 4)
      setAtomType(a, "c3");
  }
}

void GaffAtomTypeIdentifier::handleRemainingN() {
  for (auto a : atomsN_) {
    if (atomTypeSet(a))
      continue;

    verifyNeighborNumber(a, 3);

    bool hasConjugatedNeighbor = false;
    bool hasAromaticNeighbor = false;
    for (auto n : listsOfNeighbors_[a]) {
      if (isConjugated_[n])
        hasConjugatedNeighbor = true;
      if (gaffAtomTypes_[n] == "ca")
        hasAromaticNeighbor = true;
    }

    if (hasAromaticNeighbor)
      setAtomType(a, "nh");
    else if (hasConjugatedNeighbor)
      setAtomType(a, "na");
    else
      setAtomType(a, "n3");
  }
}

void GaffAtomTypeIdentifier::handleHydrogens() {
  for (auto a : atomsH_) {
    if (atomTypeSet(a))
      continue;

    verifyNeighborNumber(a, 1);

    int n = listsOfNeighbors_[a].front();
    auto nType = elementTypes_[n];

    if (nType == Utils::ElementType::P)
      setAtomType(a, "hp");
    else if (nType == Utils::ElementType::N)
      setAtomType(a, "hn");
    else if (nType == Utils::ElementType::C) {
      auto nGaffType = gaffAtomTypes_[n];
      if (nGaffType == "ca" || nGaffType == "cu" || nGaffType == "cv" || nGaffType == "c")
        setAtomType(a, "ha");
      else if (isConjugated_[n])
        setAtomType(a, "ha");
      else
        setAtomType(a, "hc");
    }
    else
      throw std::runtime_error("The atom type could not be generated for atom with index " + std::to_string(a));
  }
}

void GaffAtomTypeIdentifier::lookForCycles() {
  // Look for cycles
  for (int a = 0; a < nAtoms_; ++a)
    lookForCycle(a);

  // Check in how many cycles each atom is
  for (const auto& c : cycles_) {
    for (auto i : c) {
      ++numberOfCyclesForAtom_[i];
    }
  }
}

void GaffAtomTypeIdentifier::lookForCycle(int atomIndex) {
  /*
   * NB:
   * - comparison "XNeighbor < atomIndex" avoids duplication because cycles for
   *   lower indexes are found by lookForCycle(XNeighbor).
   * - comparison "XNeighbor > firstNeighbor" avoids duplication because otherwise each cycle is found twice.
   * - comparison "XNeighbor == (X-2)Neighbor" avoids going back to the previous atom,
   *   which was already considered in the chain.
   */
  for (auto firstNeighbor : listsOfNeighbors_[atomIndex]) {
    if (firstNeighbor < atomIndex)
      continue;

    for (auto secondNeighbor : listsOfNeighbors_[firstNeighbor]) {
      if (secondNeighbor == atomIndex || secondNeighbor < atomIndex)
        continue;

      for (auto thirdNeighbor : listsOfNeighbors_[secondNeighbor]) {
        if (thirdNeighbor == firstNeighbor || thirdNeighbor < atomIndex)
          continue;
        if (thirdNeighbor == atomIndex && secondNeighbor > firstNeighbor) {
          cycles_.push_back({atomIndex, firstNeighbor, secondNeighbor});
          continue;
        }
        for (auto fourthNeighbor : listsOfNeighbors_[thirdNeighbor]) {
          if (fourthNeighbor == secondNeighbor || fourthNeighbor < atomIndex)
            continue;
          if (fourthNeighbor == atomIndex && thirdNeighbor > firstNeighbor) {
            cycles_.push_back({atomIndex, firstNeighbor, secondNeighbor, thirdNeighbor});
            continue;
          }
          for (auto fifthNeighbor : listsOfNeighbors_[fourthNeighbor]) {
            if (fifthNeighbor == thirdNeighbor || fifthNeighbor < atomIndex)
              continue;
            if (fifthNeighbor == atomIndex && fourthNeighbor > firstNeighbor) {
              cycles_.push_back({atomIndex, firstNeighbor, secondNeighbor, thirdNeighbor, fourthNeighbor});
              continue;
            }
            for (auto sixthNeighbor : listsOfNeighbors_[fifthNeighbor]) {
              if (sixthNeighbor == fourthNeighbor || sixthNeighbor < atomIndex)
                continue;
              if (sixthNeighbor == atomIndex && fifthNeighbor > firstNeighbor) {
                cycles_.push_back({atomIndex, firstNeighbor, secondNeighbor, thirdNeighbor, fourthNeighbor, fifthNeighbor});
                continue;
              }
            }
          }
        }
      }
    }
  }
}

void GaffAtomTypeIdentifier::handleCycles() {
  for (const auto& cycle : cycles_) {
    if (cycle.size() == 6)
      handleCycle6(cycle);
    else if (cycle.size() == 5)
      handleCycle5(cycle);
    else if (cycle.size() == 4)
      handleCycle4(cycle);
    else if (cycle.size() == 3)
      handleCycle3(cycle);
  }
}

bool GaffAtomTypeIdentifier::canBeAromatic(int atomIndex) {
  // carbon
  if (elementTypes_[atomIndex] == Utils::ElementType::C && nNeighbors_[atomIndex] == 3 && gaffAtomTypes_[atomIndex] != "c")
    return true;

  // nitrogen // TODO: other possibilities?
  if (elementTypes_[atomIndex] == Utils::ElementType::N && nNeighbors_[atomIndex] == 2)
    return true;

  // TODO: sulfur, phosphorus
  return false;
}

void GaffAtomTypeIdentifier::handleCycle6(const std::list<int>& cycle) {
  // Check if aromatic:
  bool aromatic = true;
  for (auto a : cycle)
    if (!canBeAromatic(a))
      aromatic = false;

  // handle if aromatic:
  if (aromatic) {
    for (auto a : cycle) {
      if (elementTypes_[a] == Utils::ElementType::C)
        setAtomType(a, "ca", false);
      else if (elementTypes_[a] == Utils::ElementType::N)
        setAtomType(a, "nb", false);
    }
    return; // TODO: This is unnecessary, right?
  }
}

// TODO: not implemented yet
void GaffAtomTypeIdentifier::handleCycle5(const std::list<int>& /*cycle*/) {
}

void GaffAtomTypeIdentifier::handleCycle4(const std::list<int>& cycle) {
  for (auto a : cycle) {
    if (elementTypes_[a] == Utils::ElementType::C) {
      if (nNeighbors_[a] == 4)
        gaffAtomTypes_[a] = "cy";
      else if (nNeighbors_[a] == 3)
        gaffAtomTypes_[a] = "cv";
    }
  }
}

void GaffAtomTypeIdentifier::handleCycle3(const std::list<int>& cycle) {
  for (auto a : cycle) {
    if (elementTypes_[a] == Utils::ElementType::C) {
      if (nNeighbors_[a] == 4)
        gaffAtomTypes_[a] = "cx";
      else if (nNeighbors_[a] == 3)
        gaffAtomTypes_[a] = "cu";
    }
  }
}

std::string GaffAtomTypeIdentifier::conjugationReverse(const std::string& conjugatedAtomType) {
  if (conjugatedAtomType == "cc")
    return "cd";
  if (conjugatedAtomType == "ce")
    return "cf";
  if (conjugatedAtomType == "cg")
    return "ch";
  if (conjugatedAtomType == "pc")
    return "pd";
  if (conjugatedAtomType == "pe")
    return "pf";
  if (conjugatedAtomType == "nc")
    return "nd";
  if (conjugatedAtomType == "ne")
    return "nf";

  return "";
}

bool GaffAtomTypeIdentifier::atomCanConjugate(int atomIndex) {
  std::string t = gaffAtomTypes_[atomIndex];

  if (t == "c1" || t == "c2" || t == "ca" || t == "c")
    return true;

  if (t == "n2" || t == "nb" || t == "no")
    return true;

  if (t == "o")
    return true;

  if (t == "p2" || t == "p4" || t == "p5")
    return true;

  if (t == "s2" || t == "s4" || t == "s6")
    return true;

  return false;
}

std::string GaffAtomTypeIdentifier::getConjugatedType(std::string basicType, bool cyclic, bool bondAlternation) {
  if (basicType == "p4")
    return "px";
  if (basicType == "p5")
    return "py";
  if (basicType == "s4")
    return "sx";
  if (basicType == "s6")
    return "sy";

  if (basicType == "c1") {
    if (bondAlternation)
      return "cg";
    else
      return "ch";
  }

  std::string prefix;
  prefix += basicType.at(0);
  std::string suffix;
  if (cyclic) {
    if (bondAlternation)
      suffix = "c";
    else
      suffix = "d";
  }
  else {
    if (bondAlternation)
      suffix = "e";
    else
      suffix = "f";
  }
  return prefix + suffix;
}

void GaffAtomTypeIdentifier::checkConjugation() {
  std::vector<bool> canConjugate(nAtoms_, false);
  std::vector<bool> conjugationDealtWith(nAtoms_, false);

  for (int a = 0; a < nAtoms_; ++a) {
    if (atomCanConjugate(a))
      canConjugate[a] = true;
    else
      conjugationDealtWith[a] = true;
  }

  std::vector<bool> canBeTerminalForConjugation(nAtoms_, false);
  std::stack<int> terminalForConjugation;
  for (int a = 0; a < nAtoms_; ++a) {
    if (!canConjugate[a])
      continue;
    int nConjugates = 0;
    for (auto neighbor : listsOfNeighbors_[a]) {
      if (canConjugate[neighbor])
        ++nConjugates;
    }

    if (nConjugates == 1) {
      canBeTerminalForConjugation[a] = true;
      terminalForConjugation.push(a);
      isConjugated_[a] = true;
    }
  }

  while (!terminalForConjugation.empty()) {
    int startingTerminalAtom = terminalForConjugation.top();
    terminalForConjugation.pop();

    // NB: startingTerminalAtom, parent, startBond, type(c,d/e,f)
    std::stack<std::tuple<int, int, bool, bool>> toTreat;
    if (!conjugationDealtWith[startingTerminalAtom]) {
      for (auto n : listsOfNeighbors_[startingTerminalAtom]) {
        toTreat.push(std::make_tuple(n, startingTerminalAtom, false, false));
      }
      conjugationDealtWith[startingTerminalAtom] = true;
    }

    while (!toTreat.empty()) {
      auto current = toTreat.top();
      toTreat.pop();

      int newAtom = std::get<0>(current);

      // Return if it is isolated double bond
      int parent = std::get<1>(current);
      if (canBeTerminalForConjugation[parent] && canBeTerminalForConjugation[newAtom]) {
        isConjugated_[parent] = false;
        continue;
      }

      // return if already dealt with
      if (conjugationDealtWith[newAtom])
        continue;

      // Return if it is three atom with one bonded to two oxygen atoms (O=P=O) etc
      bool onlyThreeAtomsThatCanConjugate = true;
      for (auto n : listsOfNeighbors_[newAtom])
        if (canConjugate[n] && !canBeTerminalForConjugation[n])
          onlyThreeAtomsThatCanConjugate = false;
      if (onlyThreeAtomsThatCanConjugate) {
        isConjugated_[parent] = false;
        continue;
      }

      isConjugated_[newAtom] = true;

      bool bondBeginning = std::get<2>(current);
      bool bondAlternation = std::get<3>(current);
      for (auto n : listsOfNeighbors_[newAtom]) {
        toTreat.push(std::make_tuple(n, newAtom, !bondBeginning, static_cast<bool>(bondBeginning ^ bondAlternation)));
      }
      conjugationDealtWith[newAtom] = true;

      std::string currentType = gaffAtomTypes_[newAtom];
      if (currentType != "ca" && currentType != "nb" && currentType != "c" && !canBeTerminalForConjugation[newAtom]) { // NB: do not overwrite ca, nb
        std::string previousType = gaffAtomTypes_[parent];
        bool cyclic = numberOfCyclesForAtom_[newAtom] != 0;
        setAtomType(newAtom, getConjugatedType(currentType, cyclic, bondAlternation), false);
      }
    }
  }
}

void GaffAtomTypeIdentifier::readGaffAtomTypesFromFile() {
  std::ifstream indata(atomTypesFile_);
  if (!indata.is_open())
    throw std::runtime_error("The GAFF atom types file " + atomTypesFile_ + " cannot be opened.");

  std::string line;
  while (line.empty())
    std::getline(indata, line);

  int atomIndex = 0;
  while (!line.empty()) {
    std::string atomType = line;
    atomType.erase(std::remove_if(atomType.begin(), atomType.end(), ::isspace), atomType.end());
    if (atomType.length() > 2)
      throw std::runtime_error(
          "A GAFF atom type must consist either of one or two letters. Check the atom types file for errors.");
    if (atomType.length() == 0)
      break;
    std::transform(atomType.begin(), atomType.end(), atomType.begin(), ::tolower);
    if (atomIndex < nAtoms_)
      setAtomType(atomIndex, atomType);
    atomIndex++;
    if (indata.eof())
      break;
    std::getline(indata, line);
  }

  if (nAtoms_ != atomIndex)
    throw std::runtime_error("The number of atom types in the provided file does not match the number of atoms.");
}

} // namespace MolecularMechanics
} // namespace Scine
