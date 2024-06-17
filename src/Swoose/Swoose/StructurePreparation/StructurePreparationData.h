/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef PDBPREPARATION_DATAMANAGER_H
#define PDBPREPARATION_DATAMANAGER_H

#include "ProteinStructures.h"
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/NativeFilenames.h>
#include <list>
#include <vector>

namespace Scine {
namespace StructurePreparation {

/**
 * @struct StructurePreparationData StructurePreparationData.h
 * @brief This struct holds all objects used inside the structure preparation algorithm.
 */
struct StructurePreparationData {
  /**
   * @brief The structure of the system to prepare.
   */
  Utils::AtomCollection fullStructure;
  /**
   * @brief The protein substructure, i. e. all correctly identified amino acids.
   */
  Utils::AtomCollection proteinStructure;
  /**
   * @brief The nonregular substructure.
   */
  Utils::AtomCollection nonregularContainer;
  /**
   * @brief A bond order matrix of the full system.
   */
  Utils::BondOrderCollection bondOrders;
  /**
   * @brief Number of atoms in the full structure.
   */
  int numberOfAtoms{0};
  /**
   * @brief Vector of indices of the full structure that belong to the protein subsystem.
   */
  std::vector<int> vectorOfProteinIndices;
  /**
   * @brief Vector of indices of the full structure that belong to the nonRegContainer subsystem.
   */
  std::vector<int> vectorOfNonRegContainerIndices;
  /**
   * @brief The connectivity of the system. It is a vector of a list of neighbor atom indices for each atom.
   */
  std::vector<std::list<int>> listsOfNeighbors;
  /**
   * @brief Vector of all fragments containing a vector of indices that correspond to the indices
   *        of the atoms inside the fragment in the full system.
   */
  std::vector<std::vector<int>> atomIndexMapping;
  /**
   * @brief Vector containing a vector of indices that correspond to the indices
   *        of the protein and the nonregular_container in the full system.
   */
  std::vector<std::vector<int>> subsystemMapping;
  /**
   * @brief A vector of ProteinAtoms that assigns to each atom index the residue name, the atom
   * type name and the position.
   */
  std::vector<ProteinAtom> protein;
  /**
   * @brief A vector that contains the indices of all carbon atoms in the peptid bond.
   */
  std::vector<int> listOfBackboneAlphaCarbons;
};

// Default filenames
struct Default {
  std::string proteinFile = "rmc.pdb";
  std::string nonRegContainerFile = "nonregular_container.xyz";
  std::string protonatedProteinFile = "rmc_H.xyz";
  std::string protonatedNonRegContainerFile = "nonregular_container_H.xyz";
  std::string nonRegContainerInfoFile = "atomic_info_nonregular_container.dat";
  std::string systemFile = "system.xyz";
  std::string titrationSitesFile = "titrable_sites.dat";
};

struct StructurePreparationFiles {
  Default d;
  /**
   * @brief The path to the protein structure file.
   */
  std::string proteinFile;
  /**
   * @brief The path to the nonRegContainer structure file.
   */
  std::string nonRegContainerFile;
  /**
   * @brief The path to the protonated protein file. Note that per default, all pH sensitive groups are protonated
   * according to their uncharged reference state.
   */
  std::string protonatedProteinFile;
  /**
   * @brief The path to the protonated nonRegContainer file.
   */
  std::string protonatedNonRegContainerFile;
  /**
   * @brief The path to the atomic info file for the protein.
   */
  std::string atomicInfoFile;
  /**
   * @brief The path to the atomic info file for the nonRegContainer.
   */
  std::string nonRegContainerInfoFile;
  /**
   * @brief The path to which the final processed structure is written.
   */
  std::string systemFile;
  /**
   * @brief The Connectivity file that is automatically generated during the preparation.
   */
  std::string connectivityFile;
  /**
   * @brief The parent directory for all structure preparation data.
   */
  std::string preparationDataDirectory;
  /**
   * @brief The working directory.
   */
  std::string workingDirectory;
  /**
   * @brief This file contains the atom index and the residue name of a pH sensible atom.
   */
  std::string titrationSitesFile;
  void initialize() {
    // Define files here because they are needed in every mode
    proteinFile = Utils::NativeFilenames::combinePathSegments(workingDirectory, d.proteinFile);
    nonRegContainerFile = Utils::NativeFilenames::combinePathSegments(workingDirectory, d.nonRegContainerFile);
    protonatedProteinFile = Utils::NativeFilenames::combinePathSegments(workingDirectory, d.protonatedProteinFile);
    protonatedNonRegContainerFile =
        Utils::NativeFilenames::combinePathSegments(workingDirectory, d.protonatedNonRegContainerFile);
    systemFile = Utils::NativeFilenames::combinePathSegments(workingDirectory, d.systemFile);
    titrationSitesFile = Utils::NativeFilenames::combinePathSegments(workingDirectory, d.titrationSitesFile);
    nonRegContainerInfoFile = Utils::NativeFilenames::combinePathSegments(workingDirectory, d.nonRegContainerInfoFile);
  }
};

} // namespace StructurePreparation
} // namespace Scine

#endif // PDBPREPARATION_DATAMANAGER_H
