/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef STRUCTUREPREPARATION_PROTONATIONHELPER_H
#define STRUCTUREPREPARATION_PROTONATIONHELPER_H

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
struct ProtonationTypes;
namespace ProtonationHelper {

/**
 * @brief Evaluates if protonation with OpenBabel was successful.
 */
bool openBabelSuccess(std::istream& in);
/**
 * @brief Removes specific protons from an AtomCollection.
 */
void removeProtonsFromStructure(Utils::AtomCollection& structure, std::vector<int> superfluousHydrogens);
// Evaluates if an atom index is in a vector.
bool isAtomOf(const std::list<int>& group, int index);

template<std::size_t s>
bool isType(const std::array<const char*, s>& group, const std::string& atomType) {
  return (std::find(group.begin(), group.end(), atomType) != group.end());
}

} // namespace ProtonationHelper
} // namespace StructurePreparation
} // namespace Scine

#endif // STRUCTUREPREPARATION_PROTONATIONHELPER_H