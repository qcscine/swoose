/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef PDBPREPARATION_STRUCTURES_H
#define PDBPREPARATION_STRUCTURES_H

#include <Utils/Geometry/Atom.h>
#include <Utils/Typenames.h>
#include <list>
#include <map>
#include <string>

namespace Scine {
namespace Utils {
class Atom;
}
namespace StructurePreparation {

struct ProteinAtom {
  int index;
  std::string residueName;
  std::string atomType;
  Utils::Position position;
  bool isPhSensitive = false;
};

struct PeptidBond {
  int N;
  int C;
  int CA;
  int O;
};

struct ProtonationTypes {
  std::list<int> tetrahedral;
  std::list<int> pseudoTetrahedral;
  std::list<int> trigonalPlanar;
  std::list<int> linear;
};

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