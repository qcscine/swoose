/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef STRUCTUREPREPARATION_PROTONATIONHANDLER_H
#define STRUCTUREPREPARATION_PROTONATIONHANDLER_H

#include "../ProteinStructures.h"
#include "../StructurePreparationData.h"
#include "TitrationData.h"
#include <Core/Log.h>
#include <Utils/Geometry/AtomCollection.h>
#include <list>
#include <memory>

namespace Scine {
namespace Utils {
class Atom;
class Settings;
} // namespace Utils

namespace StructurePreparation {
struct StructurePreparationFiles;

class ProtonationHandler {
 public:
  /**
   * @brief Constructor
   */
  ProtonationHandler(StructurePreparationData& data, StructurePreparationFiles& files, std::shared_ptr<Utils::Settings> settings);
  ProtonationHandler();
  /**
   * @brief Protonates all amino acids
   */
  void protonateAllAminoAcids();
  /**
   * @brief Protonates the nonRegContainer using the obabel binary
   */
  void protonateNonRegContainerWithExternalSoftware();
  /**
   * @brief Changes protonation state of a pH sensitive atom.
   *
   * @return structure in which the ph sensitive amino acid is reprotonated
   */
  void reprotonateCriticalAtom(Utils::AtomCollection& structure, int indexOfCriticalAtom, const std::string& residueName);
  /*
   * @brief generate Hydrogen Atoms at the peptid chain
   */
  static Utils::Atom generatePeptidHydrogen(const Utils::AtomCollection& structure, int indexN, int indexC, int indexO);
  /**
   * @brief Set the Protein Structure object.
   */
  void setProteinStructure(Utils::AtomCollection& structure);

 private:
  void setup(const Utils::AtomCollection& structure, bool performGraphAnalysis = true);
  // Sort each atom in its correct ProtonationTypes structure
  void sortAtomTypesByProtonationCriteria();
  // Find all Peptid bonds in the structure
  void getPeptidBonds(int indexN, int& indexC, int& indexCA, int& indexO, bool& isValidPeptid);
  /*
   * @brief Protonates all sp3 hybridised groups
   */
  void protonateTetrahedralGroups(const Utils::AtomCollection& structure);
  /*
   * @brief Protonates all sp2 hybridised groups
   */
  void protonateTrigonalPlanarGroups(const Utils::AtomCollection& structure);
  /*
   * @brief Protonates all O-H and S-H groups
   */
  void protonateTerminalSAndO(const Utils::AtomCollection& structure);
  /*
   * @brief protonates the peptid bonds.
   */
  void protonatePeptidBonds(const Utils::AtomCollection& structure);
  /**
   * @brief Correctly clusters C termini into their protonation category. This can either be incomplete C termini that
   * are saturated as aldehydes or carboxylate-groups.
   */
  void detectCTerminus(const ProteinAtom& proteinAtom);

  /**
   * @brief Check whether correct OpenBabel version is present.
   */
  void checkObabelVersion();

  StructurePreparationFiles files_;
  StructurePreparationData data_;
  std::shared_ptr<Utils::Settings> settings_;
  int nAtoms_;
  ProtonationTypes protonationTypes_;
  Utils::AtomCollection hydrogenAtoms_;
  Utils::AtomCollection protein_;
  Utils::AtomCollection protonatedProtein_;

  std::vector<std::list<int>> listsOfNeighbors_;
  std::vector<int> nNeighbors_;
  std::vector<PeptidBond> peptidBonds_;

  static constexpr double minimalBondOrderToConsider_ = 0.4;
  const std::array<const char*, 5> listOfsp3AtomTypes_ = {"CB", "CG2", "CG1", "CE", "CA"};
  const std::array<const char*, 13> listOfsp2AtomTypes_ = {"CZ",  "CE1", "CE2", "CZ2", "CZ3", "CE3", "CH2",
                                                           "NE1", "NE",  "NH1", "NH2", "ND1", "ND2"};
  static constexpr std::array<const char*, 1> listOfPseudoTetrahedrals_ = {"NZ"};
  static constexpr std::array<const char*, 3> terminalOhSh_ = {"OG1", "OG", "OH"};
  static constexpr std::array<const char*, 7> AAwithsp3CG_ = {"ARG", "GLN", "GLU", "MET", "PRO", "LEU", "LYS"};
  static constexpr std::array<const char*, 4> AAwithsp3CD_ = {"ARG", "LYS", "PRO", "PYL"};
};

} // namespace StructurePreparation
} // namespace Scine

#endif // STRUCTUREPREPARATION_PROTONATIONHANDLER_H
