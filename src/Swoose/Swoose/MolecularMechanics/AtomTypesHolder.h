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
class AtomTypesHolder {
 public:
  /**
   * @brief Constructor.
   */
  explicit AtomTypesHolder(std::vector<std::string> atomTypes = {});

  /**
   * @brief Getter for the atom type for an atom with a certain index.
   */
  std::string getAtomType(int index) const;

  /**
   * @brief Returns the number of atom types stored in this object.
   */
  int size() const;

 private:
  std::vector<std::string> atomTypes_;
};

inline AtomTypesHolder::AtomTypesHolder(std::vector<std::string> atomTypes) : atomTypes_(std::move(atomTypes)) {
}

inline std::string AtomTypesHolder::getAtomType(int index) const {
  return atomTypes_.at(index);
}

inline int AtomTypesHolder::size() const {
  return atomTypes_.size();
}

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_MMATOMTYPESHOLDER_H
