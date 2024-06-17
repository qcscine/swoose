/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMMMCALCULATORSETTINGS_H
#define SWOOSE_QMMM_QMMMCALCULATORSETTINGS_H

#include <Swoose/Utilities/SettingsPopulator.h>

namespace Scine {
namespace Qmmm {

/**
 * @class QmmmCalculatorSettings QmmmCalculatorSettings.h
 * @brief Settings for QM/MM calculations.
 */
class QmmmCalculatorSettings : public Scine::Utils::Settings {
 public:
  /// @brief Populates the settings with the default settings of the given external settings object.
  void addExternalSettings(const Utils::Settings& externalSettings);
  /**
   * @brief Constructor that populates the QmmmCalculatorSettings.
   */
  QmmmCalculatorSettings() : Settings("QmmmCalculatorSettings") {
    using namespace SwooseUtilities;
    SettingsPopulator::addQmAtomsOption(_fields);
    SettingsPopulator::addElectrostaticEmbeddingOption(_fields);
    SettingsPopulator::addQmRegionXyzFileOption(_fields);
    SettingsPopulator::addIgnoreQmOption(_fields);
    SettingsPopulator::addChargeRedistributionOption(_fields);
    SettingsPopulator::addReducedQmMmEnergyOption(_fields);
    SettingsPopulator::addOptimizeLinksOption(_fields);
    SettingsPopulator::addSilenceOption(_fields);
    resetToDefaults();
  };
};

inline void QmmmCalculatorSettings::addExternalSettings(const Utils::Settings& externalSettings) {
  for (const auto& d : externalSettings.getDescriptorCollection()) {
    if (_fields.exists(d.first)) {
      continue;
    }
    _fields.push_back(d.first, d.second);
    this->addGenericValue(d.first, d.second.getDefaultValue());
  }
}

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMMMCALCULATORSETTINGS_H
