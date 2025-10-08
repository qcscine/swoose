/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_MMATOMTYPESHOLDER_H
#define MOLECULARMECHANICS_MMATOMTYPESHOLDER_H

#include <string>
#include <utility>
#include <vector>

namespace Scine {
namespace MolecularMechanics {
/**
 * @class AtomTypesHolder AtomTypesHolder.h
 * @brief Class containing the MM atom types of the atoms in a molecular system.
 */
class AtomTypesHolder : public std::vector<std::string> {
 public:
  /**
   * @brief Constructor.
   */
  explicit AtomTypesHolder(std::vector<std::string> atomTypes = {});

  /**
   * @brief Getter for the atom type for an atom with a certain index.
   */
  const std::string& getAtomType(unsigned int index) const;

  /**
   * @brief Returns a vector of unique atom types, i.e., no duplicates.
   * @return The unique atom types.
   */
  std::vector<std::string> uniqueAtomTypes() const;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_MMATOMTYPESHOLDER_H
