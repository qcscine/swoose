/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMREGIONSELECTORSETTINGS_H
#define SWOOSE_QMMM_QMREGIONSELECTORSETTINGS_H

#include <Swoose/Utilities/SettingsPopulator.h>

namespace Scine {
namespace Qmmm {

/**
 * @class QmRegionSelectorSettings QmRegionSelectorSettings.h
 * @brief Settings for automated QM region selection.
 */
class QmRegionSelectorSettings : public Utils::Settings {
 public:
  /**
   * @brief Constructor that populates the QmRegionSelectorSettings.
   */
  QmRegionSelectorSettings() : Settings("QmRegionSelectorSettings") {
    using namespace SwooseUtilities;
    SettingsPopulator::addQmRegionCenterAtom(_fields);
    SettingsPopulator::addAtomicInformationFile(_fields);
    SettingsPopulator::addParameterAndConnectivityFile(_fields, false);
    SettingsPopulator::addInitialRadius(_fields);
    SettingsPopulator::addCuttingProbability(_fields);
    SettingsPopulator::addQmRegionSizesSettings(_fields);
    SettingsPopulator::addNumAttemptsPerRadiusOption(_fields);
    SettingsPopulator::addMaxNumRefModelsOption(_fields);
    SettingsPopulator::addTolerancesForQmRegionSelection(_fields);
    SettingsPopulator::addReferenceDataGenerationOptions(_fields, false);
    SettingsPopulator::addQmRegionSelectionRandomSeed(_fields);
    SettingsPopulator::addYamlSettingsForDirectMode(_fields); // used only for internal use by app and python bindings

    // Database related settings:
    SettingsPopulator::addDatabaseSettings(_fields, "scine_swoose_qm_region_selection");
    SettingsPopulator::addDatabaseSleepTime(_fields);

    resetToDefaults();
  };
};

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMREGIONSELECTORSETTINGS_H
