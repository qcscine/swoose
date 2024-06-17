/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_STRUCTUREPREPARATIONSETTINGS_H
#define SWOOSE_STRUCTUREPREPARATIONSETTINGS_H

#include <Swoose/Utilities/SettingsPopulator.h>

namespace Scine {
namespace StructurePreparation {

/**
 * @class StructurePreparationSettings StructurePreparationSettings.h
 * @brief Settings for the PDB preparation mode.
 */
class StructurePreparationSettings : public Scine::Utils::Settings {
 public:
  /**
   * @brief Constructor that populates the StructurePreparationSettings.
   */
  StructurePreparationSettings() : Settings("StructurePreparationSettings") {
    using namespace SwooseUtilities;
    SettingsPopulator::addPreparationDataDirectory(_fields);
    SettingsPopulator::addAtomicInformationFile(_fields, false);
    SettingsPopulator::addParameterAndConnectivityFile(_fields, false);
    SettingsPopulator::addSolvation(_fields);
    SettingsPopulator::addPhValueOfSystem(_fields);
    SettingsPopulator::addChargedCandNTermini(_fields);
    SettingsPopulator::addNumberOfSolventShells(_fields);
    resetToDefaults();
  };
};

} // namespace StructurePreparation
} // namespace Scine

#endif // SWOOSE_STRUCTUREPREPARATIONSETTINGS_H
