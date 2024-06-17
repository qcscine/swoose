/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef PDBPREPARATION_PDBPREPARATIONHELPER_H
#define PDBPREPARATION_PDBPREPARATIONHELPER_H

#include "StructurePreparationData.h"
#include "boost/filesystem.hpp"
#include <Utils/Geometry/AtomCollection.h>
#include <list>
#include <memory>
#include <string>
#include <vector>

namespace Scine {
namespace Utils {
class Settings;
enum class ElementType : unsigned;
class AtomCollection;

} // namespace Utils
namespace Core {
struct Log;
}

namespace StructurePreparation {
struct ProteinAtom;
struct TitrableSite;

namespace StructurePreparationHelper {
/**
 * @brief At the beginning of each step, this function evaluates if the required files are present.
 *
 * @param mode The number of the mode that is called. 1 = prepare-analyze, 2 = prepare-protonate, 3 = prepare-finalize.
 */
void performInitialCheck(const StructurePreparationFiles& files, int mode, std::shared_ptr<Utils::Settings> settings);
/**
 * @brief This function separates the structure into a protein and a nonregular_container subsystem by performing a
 * molassembler subgraph matching with all essential amino acids.
 */
void performGraphAnalysisOnStructure(StructurePreparationData& data, const Utils::AtomCollection& structure);
/**
 * @brief Merges the substructures.
 */
void mergeProteinAndNonRegContainer(StructurePreparationData& data, const StructurePreparationFiles& files);
/**
 * @brief This function constructs each amino acid residue by recursively following the bonds along each side chain
 * and subsequently performs a subgraph matching with the amino acid that is expected.
 */
void reevaluateConnectivityForAminoAcids(StructurePreparationData& data);
/**
 * @brief After combination of (modified) substructures, the covalent bonds at the subsystem boundaries must be fixed
 * (superfluous hydrogens must be removed).
 */
void handleBoundariesBetweenProteinAndNonRegContainer(StructurePreparationData& data, StructurePreparationFiles& files);
/**
 * @brief Solvates the structure with a number of solvent shells that can be set as a setting.
 */
Utils::AtomCollection addSolvation(StructurePreparationData& data, std::shared_ptr<Utils::Settings> settings);
/**
 * @brief A method that updates several objects in the StructurePreparationData object.
 */
void updatePdbPreparationData(StructurePreparationData& data, Utils::AtomCollection& structure);
/*
 * @brief Updates the information regarding the nonRegContainer.
 */
void updateNonRegContainerVector(StructurePreparationData& data);
/**
 * @brief Searches for titrable amino acids in the structure
 */
std::vector<TitrableSite> collectTitrableSites(StructurePreparationData& data);
/**
 * @brief Removes atoms from a structure (required, if protonation state is changed).
 */
void removeAtomsFromStructure(Utils::AtomCollection& structure, std::list<int> atomsToRemove);
/**
 * @brief Evaluates, which sites are charged, and appends the corresponding atom indices to the vectors listOfPositives,
 * and listOfNegatives,
 */
void determineChargedSites(std::vector<int>& listOfNegatives, std::vector<int>& listOfPositives,
                           const StructurePreparationData& data);
/*
 * @brief This functions extracts all atoms that are part of the amino acid side chain by recursively following the
 * bonds starting from CA.
 * @param atomIndex The Index of the CA atom in the amino acid.
 * @param atomsToAdd This list is filled with all indices that belong to this amino acid.
 */
void getSideChainNeighbors(StructurePreparationData& data, int atomIndex, std::list<int>& atomsToAdd);
/*
 * @brief Detects C-Termini in the structure and moves the corresponding residue to the protein substructure.
 */
void findTermini(StructurePreparationData& data, std::list<int>& atomsToTransferToNonRegContainer);
/*
 * @brief This function maps the indices of protein and nonRegContainer to the indices in the full structure.
 */
void mapSubsystemIndicesToFullStructure(const Utils::AtomCollection& fullStructure,
                                        const Utils::AtomCollection& structure, std::vector<int>& indicesInStructure,
                                        std::vector<std::vector<int>>& subsystemMapping);
/*
 * @brief Internally transfers atoms from the protein substructure to the nonRegContainer substructure.
 */
void moveAtomsFromProteinToNonRegContainer(StructurePreparationData& data, const std::list<int>& atomsToMove);
bool compareByIndex(const ProteinAtom& a, const ProteinAtom& b);

} // namespace StructurePreparationHelper
} // namespace StructurePreparation
} // namespace Scine

#endif // PDBPREPARATION_PDBPREPARATIONHELPER_H