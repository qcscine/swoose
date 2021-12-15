/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FragmentAnalyzer.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace SwooseUtilities {

FragmentAnalyzer::FragmentAnalyzer(const std::map<int, int>& formalCharges, const std::map<int, int>& unpairedElectrons)
  : formalCharges_(formalCharges), unpairedElectrons_(unpairedElectrons) {
}

bool FragmentAnalyzer::analyzeFragment(const Utils::AtomCollection& fragment, const std::vector<int>& atomIndexMapping) {
  // Evaluate total number of electrons that should be there for a neutral fragment
  int numElectrons = 0;
  for (const auto& atom : fragment) {
    numElectrons += Utils::ElementInfo::Z(atom.getElementType());
  }

  // Evaluate charge of the fragment
  int charge = evaluateChargeOfFragment(fragment, atomIndexMapping);
  molecularCharge_ = charge;

  // Subtract the charge from the total number of electrons
  numElectrons -= charge;

  // Analyze the spin multiplicity
  return analyzeSpinMultiplicity(fragment, numElectrons, atomIndexMapping);
}

int FragmentAnalyzer::getSpinMultiplicity() const {
  return spinMultiplicity_;
}

int FragmentAnalyzer::getMolecularCharge() const {
  return molecularCharge_;
}

int FragmentAnalyzer::evaluateChargeOfFragment(const Utils::AtomCollection& fragment, const std::vector<int>& atomIndexMapping) {
  return sumUpValuesOfAtomsInFragment(fragment, atomIndexMapping, formalCharges_);
}

// TODO: Improve this function
bool FragmentAnalyzer::analyzeSpinMultiplicity(const Utils::AtomCollection& fragment, int numElectrons,
                                               const std::vector<int>& atomIndexMapping) {
  int numberUnpairedElectrons = sumUpValuesOfAtomsInFragment(fragment, atomIndexMapping, unpairedElectrons_);

  if (numberUnpairedElectrons == 0) {
    if (numElectrons % 2 == 0) {
      // No unpaired electrons expected and an even number of electrons is fine.
      spinMultiplicity_ = 1;
      return true;
    }
    else {
      // No unpaired electrons expected and an uneven number of electrons is wrong.
      return false;
    }
  }
  else if (numberUnpairedElectrons % 2 == 0) {
    if (numElectrons % 2 == 0) {
      // An even number of unpaired electrons and an even number of electrons, that is fine.
      spinMultiplicity_ = numberUnpairedElectrons + 1;
      return true;
    }
    else {
      // An even number of unpaired electrons, but an uneven number of electrons, that is wrong.
      return false;
    }
  }
  else {
    if (numElectrons % 2 == 0) {
      // If system should have an uneven number of unpaired electrons but there are an even number of electrons,
      // there is something wrong.
      return false;
    }
    else {
      // If system should have an uneven number of unpaired electrons and there are an uneven number of electrons,
      // it's fine.
      spinMultiplicity_ = numberUnpairedElectrons + 1;
      return true;
    }
  }
}

int FragmentAnalyzer::sumUpValuesOfAtomsInFragment(const Utils::AtomCollection& fragment,
                                                   const std::vector<int>& atomIndexMapping,
                                                   const std::map<int, int>& valuesMap) {
  int sum = 0;
  for (int i = 0; i < fragment.size(); ++i) {
    // Get index of atom in full system
    int indexInFullSystem = i;
    if (!atomIndexMapping.empty())
      indexInFullSystem = atomIndexMapping[i];

    // Check whether the atom is in the full system and not a saturating atom (-1 denotes that)
    if (indexInFullSystem != -1) {
      // Check whether the atom is present in the valuesMap
      if (valuesMap.find(indexInFullSystem) != valuesMap.end()) {
        sum += valuesMap.at(indexInFullSystem);
      }
    }
  }
  return sum;
}

} // namespace SwooseUtilities
} // namespace Scine
