/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ProtonationHelper.h"
#include <Utils/Constants.h>
#include <regex>

namespace Scine {
namespace StructurePreparation {
namespace ProtonationHelper {

bool openBabelSuccess(std::istream& in) {
  std::regex r1("(molecule converted)");
  std::smatch m1;

  std::string outputString;
  std::string line;
  while (std::getline(in, line))
    outputString += line;

  return std::regex_search(outputString, m1, r1);
}

void removeProtonsFromStructure(Utils::AtomCollection& structure, std::vector<int> superfluousHydrogens) {
  Utils::AtomCollection newStructure;
  for (int i = 0; i < structure.size(); ++i) {
    if (std::find(superfluousHydrogens.begin(), superfluousHydrogens.end(), i) == superfluousHydrogens.end()) {
      newStructure.push_back(structure.at(i));
    }
  }
  structure.clear();
  structure = std::move(newStructure);
}

bool isAtomOf(const std::list<int>& group, int index) {
  return std::find(group.begin(), group.end(), index) != group.end();
}

} // namespace ProtonationHelper
} // namespace StructurePreparation
} // namespace Scine