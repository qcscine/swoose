/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "StructurePreparationHelper.h"
#include "AminoAcidGraphRepresentations.h"
#include "ProteinStructures.h"
#include "Protonation/ProtonationHandler.h"
#include "Protonation/TitrationData.h"
#include "SpecialCaseHandler.h"
#include "StructurePreparationData.h"
#include "StructurePreparationSettings.h"
#include <Core/Log.h>
#include <Molassembler/Graph.h>
#include <Molassembler/Interpret.h>
#include <Molassembler/Subgraphs.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondDetector.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Solvation/SoluteSolventComplex.h>
#include <algorithm>
#include <utility>

namespace Scine {
namespace StructurePreparation {
using namespace SpecialCaseHandler;
namespace bfs = boost::filesystem;

namespace StructurePreparationHelper {

void performInitialCheck(const StructurePreparationFiles& files, int mode, std::shared_ptr<Utils::Settings> settings) {
  // Check the settings
  if (!settings->valid())
    throw Core::InitializationException("Settings are invalid");

  // Perform different initial check depending on which mode is called
  if (mode == 1) {
    return;
  }
  if (mode == 2) {
    if (!bfs::exists(files.proteinFile) && !(bfs::exists(files.nonRegContainerFile)))
      throw Core::InitializationException("Please call prepare-analyze before prepare-protonate. ");
  }
  if (mode == 3) {
    if (!bfs::exists(files.protonatedProteinFile) && !bfs::exists(files.protonatedNonRegContainerFile))
      throw Core::InitializationException("Please call prepare-protonate before prepare-finalize. ");
  }
}
std::string indexToAtomType(StructurePreparationData& data, int index) {
  for (const auto& proteinAtom : data.protein) {
    if (proteinAtom.index == index)
      return proteinAtom.atomType;
  }
}

void getSideChainNeighbors(StructurePreparationData& data, int atomIndex, std::list<int>& atomsToAdd) {
  std::string neighborName;
  auto neighbors = data.listsOfNeighbors[atomIndex];

  for (auto neighbor : neighbors) {
    neighborName = (isProteinAtom(data, neighbor)) ? indexToAtomType(data, neighbor)
                                                   : Utils::ElementInfo::symbol(data.fullStructure.getElement(neighbor));
    if (neighborName != "N" && neighborName != "C" && neighborName != "CA") {
      if (std::find(atomsToAdd.begin(), atomsToAdd.end(), neighbor) == atomsToAdd.end()) {
        atomsToAdd.push_back(neighbor);
        getSideChainNeighbors(data, neighbor, atomsToAdd);
      }
    }
    if (neighborName == "N")
      atomsToAdd.push_back(neighbor);
    else if (neighborName == "C") {
      // Save these Carbons to check for C-termini later
      data.listOfBackboneAlphaCarbons.push_back(neighbor);
      auto neighborsOfC = data.listsOfNeighbors[neighbor];
      atomsToAdd.push_back(neighbor);
      for (auto neighborOfC : neighborsOfC) {
        std::string nameOfNeighbor = (isProteinAtom(data, neighborOfC))
                                         ? indexToAtomType(data, neighborOfC)
                                         : Utils::ElementInfo::symbol(data.fullStructure.getElement(neighborOfC));
        if (nameOfNeighbor == "O" || nameOfNeighbor == "CA")
          atomsToAdd.push_back(neighborOfC);
      }
    }
  }
  // Sort vector
  atomsToAdd.sort();
  // A sanity check
  atomsToAdd.unique();
}

bool compareByIndex(const ProteinAtom& a, const ProteinAtom& b) {
  return a.index < b.index;
}

void performGraphAnalysisOnStructure(StructurePreparationData& data, const Utils::AtomCollection& structure) {
  data.vectorOfProteinIndices.clear();
  data.protein.clear();
  data.atomIndexMapping.resize(structure.size());
  Utils::BondOrderCollection bondOrders;
  bondOrders.resize(structure.size());
  bondOrders = Utils::BondDetector::detectBonds(structure);
  // perform Molassembler's graph interpretation on the full structure
  auto result = Molassembler::Interpret::graphs(structure, bondOrders);
  // Get the map of amino acid graphs
  auto mapOfAminoAcidGraphs = AminoAcids::getAminoAcidGraphs();
  // Perform graph-subgraph matching for the amino acids
  int counter = 0;
  for (long unsigned int g = 0; g < result.graphs.size(); g++) {
    counter++;
    const auto& graph = result.graphs[g];
    // iterate over the correct hierarchy of amino acids
    for (const auto* aminoAcid : AminoAcids::aminoAcidHierarchy) {
      const auto& aminoAcidGraph = mapOfAminoAcidGraphs[aminoAcid];
      // Get vectors that maps needle and haystack indices
      auto completeMappings = Molassembler::Subgraphs::complete(aminoAcidGraph, graph);
      // Iterate over all mapping vectors
      for (long unsigned int c = 0; c < completeMappings.size(); ++c) {
        for (auto iter = completeMappings[c].begin(); iter != completeMappings[c].end(); ++iter) {
          int currentIndex = data.atomIndexMapping[g][iter->right];
          // Append haystack indices to the vector of checked indices if they are new
          if (!isProteinAtom(data, currentIndex)) {
            std::string residueName = aminoAcid;
            auto residueIndex = iter->left;
            std::vector<std::string> residueTypes = AminoAcids::getResidueTypes(residueName);
            if (residueIndex > residueTypes.size())
              throw std::runtime_error("Protein and nonRegContainer separation failed");
            std::string atomType = residueTypes.at(residueIndex);
            data.vectorOfProteinIndices.push_back(currentIndex);
            // Get positions
            auto position = structure.getPosition(currentIndex) * Utils::Constants::angstrom_per_bohr;
            ProteinAtom proteinAtom;
            proteinAtom.index = currentIndex;
            proteinAtom.residueName = residueName;
            proteinAtom.atomType = atomType;
            proteinAtom.position = position;
            AminoAcidCategorizer categorizer;
            if (std::find(categorizer.critialAtomTypes.begin(), categorizer.critialAtomTypes.end(), proteinAtom.atomType) !=
                categorizer.critialAtomTypes.end())
              proteinAtom.isPhSensitive = true;
            data.protein.push_back(proteinAtom);
          }
        }
      }
    }
  }
  std::sort(data.protein.begin(), data.protein.end(), compareByIndex);
}

void updateNonRegContainerVector(StructurePreparationData& data) {
  for (int j = 0; j < data.numberOfAtoms; ++j) {
    if (!isProteinAtom(data, j)) {
      data.vectorOfNonRegContainerIndices.push_back(j);
    }
  }
}

void findTermini(StructurePreparationData& data, std::list<int>& atomsToTransferToNonRegContainer) {
  std::list<int> atomsToAdd;
  for (auto carbon : data.listOfBackboneAlphaCarbons)
    if (isCTerminus(data, carbon)) {
      auto neighbors = data.listsOfNeighbors[carbon];
      for (auto neighbor : neighbors)
        if (data.fullStructure.getElement(neighbor) == Utils::ElementType::C) {
          atomsToAdd.push_back(neighbor);
          getSideChainNeighbors(data, neighbor, atomsToAdd);
        }
    }

  for (auto atom : atomsToAdd)
    if (std::find(atomsToTransferToNonRegContainer.begin(), atomsToTransferToNonRegContainer.end(), atom) !=
        atomsToTransferToNonRegContainer.end()) {
      atomsToTransferToNonRegContainer.remove(atom);
    }
}

void moveAtomsFromProteinToNonRegContainer(StructurePreparationData& data, const std::list<int>& atomsToMove) {
  for (int index : atomsToMove) {
    data.vectorOfProteinIndices.erase(std::remove(data.vectorOfProteinIndices.begin(), data.vectorOfProteinIndices.end(), index),
                                      data.vectorOfProteinIndices.end());
    data.vectorOfNonRegContainerIndices.push_back(index);
    data.protein.erase(std::remove_if(data.protein.begin(), data.protein.end(),
                                      [&](ProteinAtom const& proteinAtom) -> bool { return proteinAtom.index == index; }),
                       data.protein.end());
  }
}

void reevaluateConnectivityForAminoAcids(StructurePreparationData& data) {
  // Initialize a list of protein indices that need to be transferred to the nonRegContainer
  // in case the mapping is zero
  std::list<int> atomsToTransferToNonRegContainer;
  for (long unsigned int i = 0; i < data.protein.size(); ++i) {
    // Start from CA to construct each residue
    if (data.protein[i].atomType == "CA") {
      std::string residueName = data.protein[i].residueName;
      Utils::AtomCollection residue;
      std::list<int> atomsToAdd;
      int initialAtom = data.protein[i].index;

      // Check connectivity by recursively following the bonds
      getSideChainNeighbors(data, initialAtom, atomsToAdd);

      for (auto a : atomsToAdd) {
        residue.push_back(data.fullStructure[a]);
      }
      // make graph comparison between residue and the expected amino acid
      auto residueBondOrders = Utils::BondDetector::detectBonds(residue);
      auto residueResult = Molassembler::Interpret::graphs(residue, residueBondOrders);
      auto residueGraph = residueResult.graphs[0];
      // Perform graph-subgraph matching
      auto mapOfAminoAcidGraphs = AminoAcids::getAminoAcidGraphs();
      for (const auto* aminoAcid : AminoAcids::aminoAcidHierarchy) {
        if (aminoAcid == residueName) {
          const auto aminoAcidGraph = mapOfAminoAcidGraphs[aminoAcid];
          // Get vectors that maps needle and haystack indices
          auto completeMappings = Molassembler::Subgraphs::complete(residueGraph, aminoAcidGraph);
          // Make atoms nonRegContainer atoms when mapping is zero
          if (completeMappings.empty()) {
            for (auto atom : atomsToAdd) {
              if (isProteinAtom(data, atom))
                atomsToTransferToNonRegContainer.push_back(atom);
            }
          }
        }
      }
    }
  }
  findTermini(data, atomsToTransferToNonRegContainer);
  atomsToTransferToNonRegContainer.sort();
  atomsToTransferToNonRegContainer.unique();
  moveAtomsFromProteinToNonRegContainer(data, atomsToTransferToNonRegContainer);
}

void updatePdbPreparationData(StructurePreparationData& data, Utils::AtomCollection& structure) {
  // Updates the StructurePreparationData object
  data.fullStructure = std::move(structure);
  data.numberOfAtoms = data.fullStructure.size();
  data.bondOrders.resize(data.numberOfAtoms);
  data.bondOrders = Utils::BondDetector::detectBonds(data.fullStructure);
  data.listsOfNeighbors =
      SwooseUtilities::TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(data.numberOfAtoms, data.bondOrders, 0.4);
  auto result = Scine::Molassembler::Interpret::graphs(data.fullStructure, data.bondOrders);
  data.atomIndexMapping.resize(data.numberOfAtoms);
  for (int i = 0; i < data.numberOfAtoms; ++i) {
    data.atomIndexMapping[result.componentMap.apply(i).component].push_back(i);
  }
}

void removeAtomsFromStructure(Utils::AtomCollection& structure, std::list<int> atomsToRemove) {
  Utils::AtomCollection newStructure;
  for (int i = 0; i < structure.size(); ++i) {
    if (std::find(atomsToRemove.begin(), atomsToRemove.end(), i) == atomsToRemove.end()) {
      newStructure.push_back(structure.at(i));
    }
  }
  structure.clear();
  structure = std::move(newStructure);
}

void mapSubsystemIndicesToFullStructure(const Utils::AtomCollection& fullStructure,
                                        const Utils::AtomCollection& structure, std::vector<int>& indicesInStructure,
                                        std::vector<std::vector<int>>& subsystemMapping) {
  for (int index = 0; index < structure.size(); ++index) {
    Utils::Atom atom = structure.at(index);
    auto indexOfAtomInFullSystem = Utils::Geometry::Distances::getIndexOfAtomInStructure(fullStructure, atom);
    indicesInStructure.push_back(indexOfAtomInFullSystem);
  }
  subsystemMapping.push_back(indicesInStructure);
}

void handleBoundariesBetweenProteinAndNonRegContainer(StructurePreparationData& data, StructurePreparationFiles& files) {
  // detect boundary regions
  std::map<int, int> boundaryAtoms;
  for (int i = 0; i < data.fullStructure.size(); ++i) {
    if (isProteinAtom(data, i) && (data.fullStructure.getElement(i) == Utils::ElementType::C)) {
      auto neighbors = data.listsOfNeighbors[i];
      for (auto n : neighbors) {
        if ((!isProteinAtom(data, n)) && (data.fullStructure.getElement(n) == Utils::ElementType::N)) {
          boundaryAtoms.insert({i, n});
        }
      }
    }
    else if (isProteinAtom(data, i) && (data.fullStructure.getElement(i) == Utils::ElementType::N)) {
      auto neighbors = data.listsOfNeighbors[i];
      for (auto n : neighbors) {
        if ((!isProteinAtom(data, n)) && (data.fullStructure.getElement(n) == Utils::ElementType::C)) {
          boundaryAtoms.insert({n, i});
        }
      }
    }
  }
  // Detect superfluous hydrogens
  std::list<int> superfluousHydrogens;
  for (auto pair : boundaryAtoms) {
    auto neighborsOfC = data.listsOfNeighbors[pair.first];
    auto neighborsOfN = data.listsOfNeighbors[pair.second];

    for (auto c : neighborsOfC)
      if (data.fullStructure.getElement(c) == Utils::ElementType::H)
        superfluousHydrogens.push_back(c);

    for (auto n : neighborsOfN)
      if (data.fullStructure.getElement(n) == Utils::ElementType::H)
        superfluousHydrogens.push_back(n);

    superfluousHydrogens.sort();
    superfluousHydrogens.unique();
  }
  // Adapt NonRegContainer Atom Collection
  std::list<int> superfluousHydrogensInNonRegContainer;
  std::list<int> superfluousHydrogensInProtein;

  for (auto hydrogen : superfluousHydrogens) {
    auto nonRegContainerIterator = std::find(data.subsystemMapping[1].begin(), data.subsystemMapping[1].end(), hydrogen);
    auto proteinIterator = std::find(data.subsystemMapping[0].begin(), data.subsystemMapping[0].end(), hydrogen);
    if (nonRegContainerIterator != data.subsystemMapping[1].end()) {
      int index = std::distance(data.subsystemMapping[1].begin(), nonRegContainerIterator);
      superfluousHydrogensInNonRegContainer.push_back(index);
    }
    if (proteinIterator != data.subsystemMapping[0].end()) {
      int index = std::distance(data.subsystemMapping[0].begin(), proteinIterator);
      superfluousHydrogensInProtein.push_back(index);
    }
  }
  removeAtomsFromStructure(data.nonregularContainer, superfluousHydrogensInNonRegContainer);
  removeAtomsFromStructure(data.proteinStructure, superfluousHydrogensInProtein);

  // Generate hydrogen atoms that need to be added
  Utils::AtomCollection hydrogenAtomsToAdd;
  for (const auto pair : boundaryAtoms) {
    int correspondingOIndex = -1;
    auto neighborsOfC = data.listsOfNeighbors[pair.first];
    for (auto neighbor : neighborsOfC) {
      if (data.fullStructure.getElement(neighbor) == Utils::ElementType::O) {
        correspondingOIndex = neighbor;
      }
    }
    if (correspondingOIndex == -1) {
      throw std::runtime_error("Generation of hydrogen atoms around a peptid bond failed.");
    }
    auto atom = Scine::StructurePreparation::ProtonationHandler::generatePeptidHydrogen(data.fullStructure, pair.second,
                                                                                        pair.first, correspondingOIndex);
    hydrogenAtomsToAdd.push_back(atom);
  }

  // remove the superfluous hydrogens
  removeAtomsFromStructure(data.fullStructure, superfluousHydrogens);
  // add the correct hydrogens
  Utils::AtomCollection totalSystem;
  totalSystem += data.proteinStructure;
  totalSystem += data.nonregularContainer;
  totalSystem += hydrogenAtomsToAdd;
  updatePdbPreparationData(data, totalSystem);

  // This final mapping is required in order to convert atomic information from the nonRegContainer to the
  // atomic info file of the entire system
  data.subsystemMapping.clear();
  data.vectorOfProteinIndices.clear();
  data.vectorOfNonRegContainerIndices.clear();
  mapSubsystemIndicesToFullStructure(data.fullStructure, data.proteinStructure, data.vectorOfProteinIndices,
                                     data.subsystemMapping);
  mapSubsystemIndicesToFullStructure(data.fullStructure, data.nonregularContainer, data.vectorOfNonRegContainerIndices,
                                     data.subsystemMapping);
}

Utils::AtomCollection addSolvation(StructurePreparationData& data, std::shared_ptr<Utils::Settings> settings) {
  // Set up the solvent
  Utils::ElementTypeCollection waterEC(3);
  waterEC.at(0) = Utils::ElementType::O;
  waterEC.at(1) = Utils::ElementType::H;
  waterEC.at(2) = Utils::ElementType::H;

  // Positions
  Utils::PositionCollection waterPC(3, 3);
  waterPC.row(0) = Utils::Position(0, 0, 0);
  waterPC.row(1) = Utils::Position(-1.63, 1.25, 0);
  waterPC.row(2) = Utils::Position(1.63, 1.25, 0);

  Utils::AtomCollection water(waterEC, waterPC);
  int numShells = settings->getInt(SwooseUtilities::SettingsNames::numberOfSolventShells);

  // The seed for choosing the surface sites
  const int seed = 5;
  Utils::SoluteSolventComplex::SolventPlacementSettings solvationSettings;
  solvationSettings.resolution = 18;
  // increment of solventOffset if no solvent molecule could be added
  solvationSettings.stepSize = 0.25;
  // number of rotamers to try for adding a solvent
  solvationSettings.numRotamers = 2;
  // maximum distance below it should be attempted to add solvent molecule
  solvationSettings.maxDistance = 10.0;
  auto solventComplex = Utils::SoluteSolventComplex::solvateShells(data.fullStructure, data.fullStructure.size(), water,
                                                                   numShells, seed, solvationSettings);
  auto mergedSolventCollection = Utils::SoluteSolventComplex::mergeSolventShellVector(solventComplex);
  return mergedSolventCollection;
}

void mergeProteinAndNonRegContainer(StructurePreparationData& data, const StructurePreparationFiles& files) {
  if (bfs::exists(files.protonatedProteinFile) && !bfs::is_empty(files.protonatedProteinFile)) {
    data.proteinStructure = Utils::ChemicalFileHandler::read(files.protonatedProteinFile).first;
  }
  if (bfs::exists(files.protonatedNonRegContainerFile) && !bfs::is_empty(files.protonatedNonRegContainerFile)) {
    data.nonregularContainer = Utils::ChemicalFileHandler::read(files.protonatedNonRegContainerFile).first;
  }
  if (data.proteinStructure.size() == 0 && data.nonregularContainer.size() == 0) {
    throw std::logic_error("Cannot parse any structure, finalization step fails!");
  }

  Utils::AtomCollection totalSystem;
  totalSystem += data.proteinStructure;
  totalSystem += data.nonregularContainer;
  updatePdbPreparationData(data, totalSystem);

  data.subsystemMapping.clear();

  mapSubsystemIndicesToFullStructure(data.fullStructure, data.proteinStructure, data.vectorOfProteinIndices,
                                     data.subsystemMapping);
  mapSubsystemIndicesToFullStructure(data.fullStructure, data.nonregularContainer, data.vectorOfNonRegContainerIndices,
                                     data.subsystemMapping);
}

std::vector<TitrableSite> collectTitrableSites(StructurePreparationData& data) {
  AminoAcidCategorizer categorizer;
  std::vector<TitrableSite> titrableSites;
  // fill data.protein object
  performGraphAnalysisOnStructure(data, data.fullStructure);
  // this code is doubled --> make this a function??
  for (long unsigned int i = 0; i < data.protein.size(); ++i) {
    // Start from CA to construct each residue
    std::string residueName = data.protein[i].residueName;
    TitrableSite site;
    if (std::find(categorizer.acids.begin(), categorizer.acids.end(), residueName) != categorizer.acids.end())
      site.isAcid = true;
    else if (std::find(categorizer.bases.begin(), categorizer.bases.end(), residueName) != categorizer.bases.end())
      site.isBase = true;

    if (data.protein[i].atomType == "CA" && (site.isAcid || site.isBase)) {
      Utils::AtomCollection residue;
      std::list<int> atomsToAdd;
      int initialAtom = data.protein[i].index;
      // Check connectivity by recursively following the bonds
      getSideChainNeighbors(data, initialAtom, atomsToAdd);
      for (auto a : atomsToAdd) {
        if (std::find(data.vectorOfNonRegContainerIndices.begin(), data.vectorOfNonRegContainerIndices.end(), a) ==
            data.vectorOfNonRegContainerIndices.end()) {
          residue.push_back(data.fullStructure[a]);
          site.indicesInFullStructure.push_back(a);
          auto it =
              std::find_if(data.protein.begin(), data.protein.end(), [a](const ProteinAtom& p) { return p.index == a; });
          if (it != data.protein.end() && it->isPhSensitive)
            site.criticalAtom = a;
        }
      }
      site.atoms = std::move(residue);
      site.residueName = residueName;
      if (site.atoms.size() > 0)
        titrableSites.push_back(site);
    }
  }
  return titrableSites;
}

void determineChargedSites(std::vector<int>& listOfNegatives, std::vector<int>& listOfPositives,
                           const StructurePreparationData& data) {
  for (int index = 0; index < data.numberOfAtoms; ++index) {
    bool negative = SpecialCaseHandler::isNegative(data, index, listOfNegatives);
    bool positive = SpecialCaseHandler::isPositive(data, index, listOfPositives);
  }
}

} // namespace StructurePreparationHelper
} // namespace StructurePreparation
} // namespace Scine
