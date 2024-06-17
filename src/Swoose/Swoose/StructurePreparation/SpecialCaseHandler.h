/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef PDBPREPARATION_SPECIALCASEHANDLER_H
#define PDBPREPARATION_SPECIALCASEHANDLER_H

#include <memory>
#include <vector>

namespace Scine {
namespace Core {
class Log;
}

namespace Utils {
class Settings;
enum class ElementType : unsigned;
} // namespace Utils

namespace StructurePreparation {
class StructurePreparationData;

namespace SpecialCaseHandler {

bool isProteinAtom(const StructurePreparationData& data, int index);
/**
 * @brief Evaluates whether an atom with given index bears a negative charge.
 *
 * @param data The StructurePreparationData object given as a reference.
 * @param index The index of the checked atom.
 * @param listOfNegatives A vector of negatively charged atoms updated iteratively.
 * @return true if atom with index has a charge of -1
 */
bool isNegative(const StructurePreparationData& data, int index, std::vector<int>& listOfNegatives);
/**
 * @brief Evaluates whether an atom with given index bears a positive charge.
 *
 * @param data The StructurePreparationData object given as a reference.
 * @param index The index of the atom for which the positive charge is checked.
 * @param listOfPositives A vector of positively charged atoms updated iteratively.
 * @return true if atom with index has a charge of +1
 */
bool isPositive(const StructurePreparationData& data, int index, std::vector<int>& listOfPositives);
/**
 * @brief Evaluates if the carbon atom with given index is a C-Terminus
 * @return true If the atom with the corresponding index is both a protein atom and a carboxylate Atom bounded to CA.
 */
bool isCTerminus(const StructurePreparationData& data, int index);
/**
 * @brief Evaluates if the carbon atom with given index is a Carboxylate-C.
 */
bool isCarboxylateC(const StructurePreparationData& data, int index);

} // namespace SpecialCaseHandler
} // namespace StructurePreparation
} // namespace Scine

#endif // PDBPREPARATION_SPECIALCASEHANDLER_H