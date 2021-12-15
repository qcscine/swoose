/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_SFAMCALCULATORSETTINGS_H
#define MOLECULARMECHANICS_SFAMCALCULATORSETTINGS_H

#include <Swoose/Utilities/SettingsPopulator.h>

namespace Scine {
namespace MolecularMechanics {

/**
 * @class SfamCalculatorSettings SfamCalculatorSettings.h
 * @brief Settings for SFAM molecular mechanics calculations.
 */
class SfamCalculatorSettings : public Scine::Utils::Settings {
 public:
  /**
   * @brief Constructor that populates the SfamCalculatorSettings.
   */
  SfamCalculatorSettings() : Settings("SfamCalculatorSettings") {
    using namespace SwooseUtilities;
    SettingsPopulator::addSfamAtomTypeLevel(_fields);
    SettingsPopulator::addPrintContributionsMolecularMechanicsOption(_fields);
    SettingsPopulator::addOnlyCalculateBondedContribution(_fields);
    SettingsPopulator::addParameterAndConnectivityFile(_fields);
    SettingsPopulator::addDetectBondsWithCovalentRadiiOption(_fields);
    SettingsPopulator::addNonCovalentCutoffRadius(_fields);
    SettingsPopulator::addHydrogenBondCorrection(_fields);
    SettingsPopulator::addApplyCutoffDuringInitializationOption(_fields);
    resetToDefaults();
  };
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_SFAMCALCULATORSETTINGS_H
