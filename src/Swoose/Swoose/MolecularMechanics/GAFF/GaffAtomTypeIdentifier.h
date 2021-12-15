/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_GAFFATOMTYPEIDENTIFIER_H
#define MOLECULARMECHANICS_GAFFATOMTYPEIDENTIFIER_H

#include "../AtomTypesHolder.h"
#include <Utils/Typenames.h>
#include <list>
#include <string>
#include <vector>

namespace Scine {
namespace MolecularMechanics {

/**
 * @class GaffAtomTypeIdentifier GaffAtomTypeIdentifier.h
 * @brief Class determining the GAFF atom types of a system given its connectivity,
 *        its element types and an atom type level.
 */
class GaffAtomTypeIdentifier {
 public:
  /**
   * @brief Constructor.
   * @param nAtoms Number of atoms.
   * @param elementTypes The elements.
   * @param listsOfNeighbors The connectivity.
   * @param atomTypesFile File that holds the GAFF atom types. If empty (as default), use automated algorithm instead.
   */
  GaffAtomTypeIdentifier(int nAtoms, Utils::ElementTypeCollection elementTypes,
                         const std::vector<std::list<int>>& listsOfNeighbors, std::string atomTypesFile = "");

  /**
   * @brief Set the atom type for an atom.
   * @param index The atom's index
   * @param type The atom type to set for that atom.
   * @param throwIfAlreadySet If this parameter is true and the atom type was already set, an exception is thrown.
   */
  void setAtomType(int index, std::string type, bool throwIfAlreadySet = true);
  /**
   * @brief This function determines the GAFF atom types and returns an instance of AtomTypeHolder.
   */
  AtomTypesHolder getAtomTypes();

 private:
  bool atomTypeSet(int atomIndex) const;
  void setArraysForElementTypes();
  void verifyNeighborNumber(int atomIndex, int imposedNumberOfNeighbors);
  void handleHalogens();
  void handleOxygenFunctionalGroups();
  void handleSulfurFunctionalGroups();
  void handleSpecialNitrogens();
  void handlePhosphorus();
  void handleCarbons();
  void handleRemainingN();
  void handleHydrogens();
  void calculateNumberOfNeighbors();
  void lookForCycles();
  void lookForCycle(int atomIndex);
  void handleCycles();
  bool canBeAromatic(int atomIndex);
  void handleCycle6(const std::list<int>& cycle);
  void handleCycle5(const std::list<int>& cycle);
  void handleCycle4(const std::list<int>& cycle);
  void handleCycle3(const std::list<int>& cycle);
  std::string conjugationReverse(const std::string& conjugatedAtomType);
  bool atomCanConjugate(int atomIndex);
  std::string getConjugatedType(std::string basicType, bool cyclic, bool bondAlternation);
  void checkConjugation();

  int nAtoms_;
  Utils::ElementTypeCollection elementTypes_;
  const std::vector<std::list<int>>& listsOfNeighbors_;
  std::vector<int> nNeighbors_;
  std::vector<int> numberOfCyclesForAtom_;
  std::list<std::list<int>> cycles_;
  std::vector<std::string> gaffAtomTypes_;
  std::vector<bool> isConjugated_;
  std::vector<int> atomsH_, atomsC_, atomsN_, atomsO_, atomsP_, atomsS_, atomsF_, atomsCl_, atomsBr_, atomsI_;

  /*
   * Function to read the atom types from file.
   */
  void readGaffAtomTypesFromFile();
  /*
   * File which holds the GAFF atom types for the system (manual input).
   * If empty string, use the normal algorithm to determine
   */
  std::string atomTypesFile_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_GAFFATOMTYPEIDENTIFIER_H
