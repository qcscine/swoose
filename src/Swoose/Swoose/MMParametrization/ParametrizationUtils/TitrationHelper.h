/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef PDBPREPARATION_TITRATIONHELPER_H
#define PDBPREPARATION_TITRATIONHELPER_H

#include <memory>
#include <string>
#include <vector>

namespace Scine {
namespace Utils {
class Atom;
class Settings;
class AtomCollection;
} // namespace Utils

namespace Core {
struct Log;
}

namespace StructurePreparation {
struct TitrableSite;
} // namespace StructurePreparation

namespace MMParametrization {
struct TrainingData;
struct TitrationResults;

class TitrationHelper {
 public:
  TitrationHelper(std::shared_ptr<Utils::Settings>& settings);
  /**
   * @brief Changes the protonation state of a pH sensitive site in a fragment during the parametrization.
   *
   * @param refStructure The input structure.
   * @param residueName The residue name.
   * @param isBase Whether the amino acid is a base.
   * @param indexOfCriticalAtom The index of the atom at which the protonation pattern should be changed.
   * @param superfluousHydrogens The hydrogens that must be eliminated.
   * @return Utils::AtomCollection The updated structure.
   */
  Utils::AtomCollection changeProtonationState(const Utils::AtomCollection& refStructure,
                                               const std::string& residueName, bool isBase, int indexOfCriticalAtom,
                                               const std::vector<int>& superfluousHydrogens);

  /**
   * @brief Read training data from a directory.
   * @param results The TitrationResults object defining which functional groups are expected in the training data
   *                directory and the path to the directory. Furthermore, the training data is saved in the results
   *                object.
   */
  void collectTrainingData(TitrationResults& results);
  /**
   * @brief Estimate the free energy for each protonation site. This function updates the data in the results object.
   * @param results The results object defining the protonation sites.
   */
  static void calculateFreeEnergiesOfDeprotonation(TitrationResults& results);
  /**
   * @brief Try to determine the pKa value of the amino-acid residue with the given name. An exception is thrown if the
   *        residue is not tabulated,
   * @param residueName the name of the residue.
   * @return The pKa value.
   */
  double getModelPka(const std::string& residueName);

 private:
  double getPkaOfTrainingMolecule(std::string dataDirectory);
  double getEnergyOfDeprotonationForTrainingMolecule(std::string dataDirectory);

  // The settings
  std::shared_ptr<Utils::Settings>& settings_;

  const std::vector<std::string> availableFunctionalGroups_ = {"Phenol", "Alcohol", "NH3", "SH", "Imidazole", "COOH"};
};
} // namespace MMParametrization
} // namespace Scine

#endif // PDBPREPARATION_TITRATIONHELPER_H
