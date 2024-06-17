/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_MMPARAMETRIZATIONSETTINGS_H
#define MMPARAMETRIZATION_MMPARAMETRIZATIONSETTINGS_H

#include <Swoose/Utilities/SettingsPopulator.h>

namespace Scine {
namespace MMParametrization {

/**
 * @class MMParametrizationSettings MMParametrizationSettings.h
 * @brief Settings for the MM model parametrizer.
 */
class MMParametrizationSettings : public Scine::Utils::Settings {
 public:
  /**
   * @brief Constructor that populates the MMParametrizationSettings.
   */
  MMParametrizationSettings() : Settings("MMParametrizationSettings") {
    using namespace SwooseUtilities;
    SettingsPopulator::addSfamAtomTypeLevel(_fields);
    SettingsPopulator::addBondOrderThreshold(_fields);
    SettingsPopulator::addConnectivityRefinementOption(_fields);
    SettingsPopulator::addParameterAndConnectivityFile(_fields);
    SettingsPopulator::addExistingParameterFile(_fields);
    SettingsPopulator::addConstrainMMParametersOption(_fields);
    SettingsPopulator::addOptimizeImproperDihedralForceConstantsOption(_fields);
    SettingsPopulator::addNumberAtomsThreshold(_fields);
    SettingsPopulator::addSubsystemRadius(_fields);
    SettingsPopulator::addReferenceDataGenerationOptions(_fields);
    SettingsPopulator::addDatabaseSettings(_fields, "scine_swoose_mm_parametrization");
    SettingsPopulator::addDatabaseSleepTime(_fields);
    SettingsPopulator::addAtomicInformationFile(_fields);
    SettingsPopulator::addEarlyTerminationOption(_fields);
    SettingsPopulator::addReuseDatabaseOption(_fields);
    SettingsPopulator::addTerminateAfterReferenceDataGenerationOption(_fields);
    SettingsPopulator::addUseCsvInputFormatOption(_fields);
    SettingsPopulator::addConvertToCm5Option(_fields);
    SettingsPopulator::addTitration(_fields);
    SettingsPopulator::addUseThermochemistryForTitration(_fields);
    SettingsPopulator::addTrainingDataDirectory(_fields);
    SettingsPopulator::addYamlSettingsForDirectMode(_fields); // used only for internal use by app and python bindings
    SettingsPopulator::addTitrationSiteFile(_fields);

    // For Gaussian, Orca and Turbomole calculations
    SettingsPopulator::addUseGaussianOption(_fields);
    SettingsPopulator::addBaseWorkingDirectory(_fields);
    SettingsPopulator::addGaussianMethodAndBasisSet(_fields);
    SettingsPopulator::addExternalProgramNProcs(_fields);

    // For QM calculators in general
    SettingsPopulator::addReferenceProgram(_fields);
    SettingsPopulator::addReferenceMethodAndBasisSet(_fields);
    SettingsPopulator::addIncreaseScfSafetyOption(_fields);

    resetToDefaults();
  };
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_MMPARAMETRIZATIONSETTINGS_H
