/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_SFAMATOMTYPEIDENTIFIER_H
#define MOLECULARMECHANICS_SFAMATOMTYPEIDENTIFIER_H

#include "../AtomTypesHolder.h"
#include <Utils/Typenames.h>
#include <list>
#include <string>
#include <vector>

namespace Scine {
namespace MolecularMechanics {

/// @brief Atom type level
enum class SfamAtomTypeLevel { Elements, Low, High, Unique };

/**
 * @class AtomTypeIdentifier AtomTypeIdentifier.h
 * @brief Class determining the atom types of a system given its connectivity, its element types and an atom type level.
 */
class SfamAtomTypeIdentifier {
 public:
  /**
   * @brief Constructor.
   * @param nAtoms Number of atoms.
   * @param elementTypes The elements.
   * @param listsOfNeighbors The connectivity.
   */
  SfamAtomTypeIdentifier(int nAtoms, Utils::ElementTypeCollection elementTypes,
                         const std::vector<std::list<int>>& listsOfNeighbors);

  /**
   * @brief Set the atom type for an atom.
   */
  void setAtomType(int index, std::string type);
  /**
   * @brief This function determines the atom types according to the atom type level and returns an instance
   *        of AtomTypeHolder.
   */
  AtomTypesHolder getAtomTypes(SfamAtomTypeLevel atl);
  /**
   * @brief This function generates an atom type level from its string representation.
   */
  static SfamAtomTypeLevel generateSfamAtomTypeLevelFromString(std::string sfamAtomTypeLevelString);
  /**
   * @brief Read the SFAM atom types from file.
   * @param sfamAtomTypesFileName Name of the file containing the SFAM atom types.
   * @param nAtoms The number of atoms in the system. This value is used for sanity checking that every atom has its
   * type.
   * @return  The atom types.
   */
  static AtomTypesHolder getAtomTypesFromFile(const std::string& sfamAtomTypesFileName, unsigned int nAtoms);
  /**
   * @brief Write the SFAM atom types to a file.
   * @param fileName The name of the file to which the atom types are written.
   * @param atomTypes The atom types.
   */
  static void writeAtomTypesToFile(const std::string& fileName, const AtomTypesHolder& atomTypes);

 private:
  int nAtoms_;
  Utils::ElementTypeCollection elementTypes_;
  const std::vector<std::list<int>>& listsOfNeighbors_;
  std::vector<std::string> atomTypes_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_SFAMATOMTYPEIDENTIFIER_H
