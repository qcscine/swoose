/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef PDBPREPARATION_PDBPROCESSOR_H
#define PDBPREPARATION_PDBPROCESSOR_H

#include "StructurePreparationData.h"
#include <Core/Log.h>
#include <string>

namespace Scine {

namespace Core {
struct Log;
}
namespace Utils {
class Settings;
}

namespace StructurePreparation {
class ProtonationHandler;

class StructureProcessor {
 public:
  static constexpr const char* model = "PDB_Processor";
  /**
   * @brief Constructor.
   */
  StructureProcessor();
  /**
   * @brief Main function of this class. It prepares the protein input structure.
   * @param structure The input file.
   * @param mode The mode. Valid modes are: prepare-analyze, prepare-protonate and prepare-finalize.
   */
  void prepare(const std::string& structureFile, const std::string& mode);
  /**
   * @brief Accessor for the settings.
   */
  Utils::Settings& settings();
  /**
   * @brief Constant accessor for the settings.
   */
  const Utils::Settings& settings() const;
  /**
   * @brief Getter for the name of the Parametrizer.
   */
  std::string name() const;
  /**
   * @brief Analyzes the input structure.
   */
  void analyzeStructure(const std::string& structureFile);
  /**
   * @brief Protonates protein and nonRegContainer separately.
   */
  void protonate();
  /**
   * @brief This function merges protein and nonRegContainer, corrects the boundaries, adds optional solvation and
   * generates atomic info file and connectivity file.
   */
  void finalize();
  // Sets the logger.
  void setLog(Core::Log& log);
  // Sets the corresponding files.
  void setFiles(StructurePreparationFiles files);

 private:
  /*
   * @brief Reads in the structure from either a PDB file or an XYZ file.
   *
   * @param includeH Also parse hydrogen atoms.
   * @return Utils::AtomCollection The AtomCollection of the parsed structure.
   */
  Utils::AtomCollection getAtomCollectionFromInput(const std::string& structureFile, bool includeH) const;
  /*
   * @brief This function sets the paths to the working directory and all files that are needed in each mode.
   */
  void performInitialSetup(const std::string& structureFile);
  /*
   * @brief Checks if nonRegContainer is present. If no, directly proceeds to protonation and finalization. If yes,
   * writes the nonRegContainer substructure to an XYZ file and stops.
   */
  void processNonRegContainerSubstructure();
  // The logger.
  Core::Log log_;
  // This object holds all the data used in the structure preparation algorithm.
  StructurePreparationData data_;
  // This object holds all the files needed for the structure Preparation.
  StructurePreparationFiles files_;
  // The settings
  std::shared_ptr<Utils::Settings> settings_;
  //
  bool protonationDone_ = false;
  bool finalizationDone_ = false;
};

} // namespace StructurePreparation
} // namespace Scine

#endif // PDBPREPARATION_PdbProcessor_H
