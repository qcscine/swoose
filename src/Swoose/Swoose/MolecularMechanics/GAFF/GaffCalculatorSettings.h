/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_GAFFCALCULATORSETTINGS_H
#define MOLECULARMECHANICS_GAFFCALCULATORSETTINGS_H

#include <Swoose/Utilities/SettingsPopulator.h>

namespace Scine {
namespace MolecularMechanics {

/**
 * @class GaffCalculatorSettings GaffCalculatorSettings.h
 * @brief Settings for GAFF molecular mechanics calculations.
 */
class GaffCalculatorSettings : public Scine::Utils::Settings {
 public:
  /**
   * @brief Constructor that populates the GaffCalculatorSettings.
   */
  GaffCalculatorSettings() : Settings("GaffCalculatorSettings") {
    using namespace SwooseUtilities;
    SettingsPopulator::addPrintContributionsMolecularMechanicsOption(_fields);
    SettingsPopulator::addOnlyCalculateBondedContribution(_fields);
    SettingsPopulator::addParameterAndConnectivityFile(_fields, true);
    SettingsPopulator::addDetectBondsWithCovalentRadiiOption(_fields);
    SettingsPopulator::addNonCovalentCutoffRadius(_fields);
    SettingsPopulator::addGaffAtomicChargesFile(_fields);
    SettingsPopulator::addGaffAtomTypesFile(_fields);
    SettingsPopulator::addApplyCutoffDuringInitializationOption(_fields);
    resetToDefaults();
  };
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_GAFFCALCULATORSETTINGS_H
