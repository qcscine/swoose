/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Swoose/MolecularMechanics/AtomTypesHolder.h"
#include <set>

namespace Scine {
namespace MolecularMechanics {

AtomTypesHolder::AtomTypesHolder(std::vector<std::string> atomTypes) : std::vector<std::string>(std::move(atomTypes)) {
}

const std::string& AtomTypesHolder::getAtomType(unsigned int index) const {
  return this->at(index);
}

std::vector<std::string> AtomTypesHolder::uniqueAtomTypes() const {
  std::set<std::string> uniqueAtomTypes;
  for (const auto& atomType : *this) {
    uniqueAtomTypes.insert(atomType);
  }
  std::vector<std::string> uniqueAtomTypeVector;
  for (const auto& atomType : uniqueAtomTypes) {
    uniqueAtomTypeVector.push_back(atomType);
  }
  return uniqueAtomTypeVector;
}

} // namespace MolecularMechanics
} // namespace Scine
