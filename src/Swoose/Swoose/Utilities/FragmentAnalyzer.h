/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_FRAGMENTANALYZER_H
#define SWOOSEUTILITIES_FRAGMENTANALYZER_H

#include <map>
#include <vector>

namespace Scine {

namespace Utils {
class AtomCollection;
} // namespace Utils

namespace SwooseUtilities {

class FragmentAnalyzer {
 public:
  /**
   * @brief Constructor.
   * @param formalCharges A map containing the indices of atoms and their formal charge in the full system.
   * @param unpairedElectrons A map containing the indices of atoms and the number of unpaired electrons that can
   *                          be assigned to that atom.
   */
  FragmentAnalyzer(const std::map<int, int>& formalCharges, const std::map<int, int>& unpairedElectrons);
  /**
   * @brief Determines whether a given fragment is valid and if it is,
   *        evaluates the charge and multiplicity of the fragment.
   * @param fragment The molecular structure representing the fragment.
   * @param atomIndexMapping Vector of indices that correspond to the indices of the atoms
   *                         inside the given fragment in the full system. Default value is an empty vector,
   *                         i.e. the indices of the fragment equal the full system's indices.
   * @return bool Returns whether fragment is valid.
   */
  bool analyzeFragment(const Utils::AtomCollection& fragment, const std::vector<int>& atomIndexMapping = {});
  /**
   * @brief Getter for the molecular charge of the fragment that was previously analyzed.
   */
  int getMolecularCharge() const;
  /**
   * @brief Getter for the spin multiplicity of the fragment that was previously analyzed.
   */
  int getSpinMultiplicity() const;

 private:
  /*
   * @brief Evaluates the charge of a fragment based on the formal charges map in the ParametrizationData object.
   */
  int evaluateChargeOfFragment(const Utils::AtomCollection& fragment, const std::vector<int>& atomIndexMapping);
  /*
   * @brief Analyzes the spin multiplicity, updates the corresponding member variable
   *        and returns whether the fragment is valid.
   */
  bool analyzeSpinMultiplicity(const Utils::AtomCollection& fragment, int numElectrons, const std::vector<int>& atomIndexMapping);
  /*
   * @brief Sums up either the formal charges or the unpaired electrons of all atoms in a given fragment.
   *        The information about which atoms have which values is given through a map.
   */
  int sumUpValuesOfAtomsInFragment(const Utils::AtomCollection& fragment, const std::vector<int>& atomIndexMapping,
                                   const std::map<int, int>& valuesMap);
  // The molecular charge of the fragment that was previously analyzed.
  int molecularCharge_ = 0;
  // The spin multiplicity of the fragment that was previously analyzed.
  int spinMultiplicity_ = 1;
  /*
   * @brief A map containing the indices of atoms and their formal charge in the full system.
   */
  const std::map<int, int>& formalCharges_;
  /*
   * @brief A map containing the indices of atoms and the number of unpaired electrons that can
   *        be assigned to that atom.
   */
  const std::map<int, int>& unpairedElectrons_;
};

} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_FRAGMENTANALYZER_H
