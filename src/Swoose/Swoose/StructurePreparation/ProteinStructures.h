/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef PDBPREPARATION_STRUCTURES_H
#define PDBPREPARATION_STRUCTURES_H

#include <Utils/Typenames.h>
#include <list>
#include <map>
#include <string>

namespace Scine {
namespace Utils {
class Atom;
}
namespace StructurePreparation {
/**
 * @brief This struct represents an atom in a protein, characterized by its index, residue name, atom type, etc.
 */
struct ProteinAtom {
  int index = std::numeric_limits<int>::infinity();
  std::string residueName;
  std::string atomType;
  Utils::Position position;
  bool isPhSensitive = false;
};

/**
 * @brief This struct collects the atom indices involved in a peptide bond.
 */
struct PeptidBond {
  int N = std::numeric_limits<int>::infinity();
  int C = std::numeric_limits<int>::infinity();
  int CA = std::numeric_limits<int>::infinity();
  int O = std::numeric_limits<int>::infinity();
};

struct ProtonationTypes {
  std::list<int> tetrahedral;
  std::list<int> pseudoTetrahedral;
  std::list<int> trigonalPlanar;
  std::list<int> linear;
};

/**
 * This struct collects data for some amino acids with their pKa value.
 */
struct AminoAcidCategorizer {
  std::vector<std::string> acids = {"ASP", "GLU", "CYS", "TYR"}; // in ref state protonated
  std::vector<std::string> bases = {"ARG", "HIS", "LYS"};        // in ref state deprotonated
  std::vector<std::string> critialAtomTypes = {"OD2", "OE2", "SG", "OH", "NH1", "NE2", "NZ"};
  std::map<std::string, double> modelPkaMap = {{"ASP", 4.0},  {"GLU", 4.4}, {"CYS", 9.5}, {"TYR", 9.6},
                                               {"ARG", 12.0}, {"HIS", 7.0}, {"LYS", 10.4}};
  std::map<std::string, std::string> functionalGroups = {{"ASP", "COOH"},   {"GLU", "COOH"}, {"CYS", "SH"},
                                                         {"TYR", "Phenol"}, {"ARG", "NH3"},  {"HIS", "Imidazole"},
                                                         {"LYS", "NH3"}};
};

} // namespace StructurePreparation
} // namespace Scine

#endif // PDBPREPARATION_STRUCTURES_H
