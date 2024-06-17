/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_SETTINGSPOPULATOR_H
#define SWOOSEUTILITIES_SETTINGSPOPULATOR_H

#include "OptionNames.h"
#include "SettingsNames.h"
#include <Utils/ExternalQC/Gaussian/GaussianCalculatorSettings.h>
#include <Utils/ExternalQC/Orca/OrcaCalculatorSettings.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace SwooseUtilities {

/**
 * @class SettingsPopulator SettingsPopulator.h
 * @brief Class with only static functions to populate settings.
 */
class SettingsPopulator {
 public:
  // These functions used mostly for MM
  static void addSfamAtomTypeLevel(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addPrintContributionsMolecularMechanicsOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addOnlyCalculateBondedContribution(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addParameterAndConnectivityFile(Utils::UniversalSettings::DescriptorCollection& settings,
                                              bool setEmptyDefault = true);
  static void addDetectBondsWithCovalentRadiiOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addNonCovalentCutoffRadius(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addHydrogenBondCorrection(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addApplyCutoffDuringInitializationOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addGaffAtomicChargesFile(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addGaffAtomTypesFile(Utils::UniversalSettings::DescriptorCollection& settings);

  // These functions used for QM/MM
  static void addQmAtomsOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addElectrostaticEmbeddingOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addQmRegionXyzFileOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addIgnoreQmOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addOptimizeLinksOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addSilenceOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addChargeRedistributionOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addReducedQmMmEnergyOption(Utils::UniversalSettings::DescriptorCollection& settings);

  // These functions are mostly for the MM parametrization
  static void addBondOrderThreshold(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addConnectivityRefinementOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addExistingParameterFile(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addConstrainMMParametersOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addOptimizeImproperDihedralForceConstantsOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addGaussianMethodAndBasisSet(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addExternalProgramNProcs(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addBaseWorkingDirectory(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addNumberAtomsThreshold(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addSubsystemRadius(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addReferenceDataGenerationOptions(Utils::UniversalSettings::DescriptorCollection& settings,
                                                bool includeReadWriteMode = true);
  static void addDatabaseSettings(Utils::UniversalSettings::DescriptorCollection& settings, std::string defaultDatabaseName);
  static void addDatabaseSleepTime(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addAtomicInformationFile(Utils::UniversalSettings::DescriptorCollection& settings, bool setEmptyDefault = true);
  static void addUseGaussianOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addReferenceProgram(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addReferenceMethodAndBasisSet(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addIncreaseScfSafetyOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addEarlyTerminationOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addReuseDatabaseOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addTerminateAfterReferenceDataGenerationOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addUseCsvInputFormatOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addConvertToCm5Option(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addYamlSettingsForDirectMode(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addTitration(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addUseThermochemistryForTitration(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addTrainingDataDirectory(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addTitrationSiteFile(Utils::UniversalSettings::DescriptorCollection& settings);

  // These functions used for QM region selection
  static void addQmRegionCenterAtoms(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addInitialRadius(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addCuttingProbability(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addQmRegionSizesSettings(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addNumAttemptsPerRadiusOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addMaxNumRefModelsOption(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addTolerancesForQmRegionSelection(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addQmRegionSelectionRandomSeed(Utils::UniversalSettings::DescriptorCollection& settings);

  // These functions used for QM region selection via a Database
  static void addMethodFamily(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addProgram(Utils::UniversalSettings::DescriptorCollection& settings);

  // These functions used for structure preparation
  static void addPreparationDataDirectory(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addSolvation(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addPhValueOfSystem(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addNumberOfSolventShells(Utils::UniversalSettings::DescriptorCollection& settings);
  static void addChargedCandNTermini(Utils::UniversalSettings::DescriptorCollection& settings);
  /**
   * @brief Deleted constructor.
   */
  SettingsPopulator() = delete;
};

inline void SettingsPopulator::addSfamAtomTypeLevel(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor sfamAtomTypeLevel("Sets the atom type level for SFAM's MM model.");
  sfamAtomTypeLevel.addOption("elements");
  sfamAtomTypeLevel.addOption("low");
  sfamAtomTypeLevel.addOption("high");
  sfamAtomTypeLevel.addOption("unique");
  sfamAtomTypeLevel.setDefaultOption("high");
  settings.push_back(SettingsNames::sfamAtomTypeLevel, std::move(sfamAtomTypeLevel));
}

inline void
SettingsPopulator::addPrintContributionsMolecularMechanicsOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor printContributionsMolecularMechanics(
      "Sets the option to have a very verbose output from the MM calculation, which includes the individual energy "
      "contributions.");
  printContributionsMolecularMechanics.setDefaultValue(false);
  settings.push_back(SettingsNames::printContributionsMolecularMechanics, std::move(printContributionsMolecularMechanics));
}

inline void SettingsPopulator::addOnlyCalculateBondedContribution(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor onlyCalculateBondedContribution(
      "Sets the option to only calculate covalent contributions within the MM model.");
  onlyCalculateBondedContribution.setDefaultValue(false);
  settings.push_back(SettingsNames::onlyCalculateBondedContribution, std::move(onlyCalculateBondedContribution));
}

inline void SettingsPopulator::addParameterAndConnectivityFile(Utils::UniversalSettings::DescriptorCollection& settings,
                                                               bool setEmptyDefault) {
  Utils::UniversalSettings::StringDescriptor parameterFilePath(
      "Path to the MM parameter file (for reading and writing).");
  if (setEmptyDefault)
    parameterFilePath.setDefaultValue("");
  else
    parameterFilePath.setDefaultValue("Parameters.dat");
  settings.push_back(Utils::SettingsNames::parameterFilePath, std::move(parameterFilePath));

  Utils::UniversalSettings::StringDescriptor connectivityFilePath(
      "Path to the system connectivity file (for reading and writing).");
  if (setEmptyDefault)
    connectivityFilePath.setDefaultValue("");
  else
    connectivityFilePath.setDefaultValue("Connectivity.dat");
  settings.push_back(SettingsNames::connectivityFilePath, std::move(connectivityFilePath));
}

inline void SettingsPopulator::addDetectBondsWithCovalentRadiiOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor detectBonds(
      "Decides whether the connectivity should be determined by bond detection based on covalent radii instead of "
      "reading the connectivity file.");
  detectBonds.setDefaultValue(false);
  settings.push_back(SettingsNames::detectBondsWithCovalentRadii, std::move(detectBonds));
}

inline void SettingsPopulator::addNonCovalentCutoffRadius(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor cutoffRadius(
      "The cutoff radius for non covalent interactions in Angstrom.");
  cutoffRadius.setMinimum(0.0);
  cutoffRadius.setDefaultValue(1200.0);
  settings.push_back(SettingsNames::nonCovalentCutoffRadius, std::move(cutoffRadius));
}

inline void SettingsPopulator::addHydrogenBondCorrection(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor hydrogenBondCorrection("Include hydrogen bond interaction in MM model.");
  hydrogenBondCorrection.setDefaultValue(true);
  settings.push_back(SettingsNames::hydrogenBondCorrection, std::move(hydrogenBondCorrection));
}

inline void SettingsPopulator::addApplyCutoffDuringInitializationOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor applyCutoffDuringInitialization(
      "Decides whether the non-covalent cutoff radius should be enforced during the initialization of the calculator "
      "to "
      "exclude the interactions beyond the distance cutoff.");
  applyCutoffDuringInitialization.setDefaultValue(false);
  settings.push_back(SettingsNames::applyCutoffDuringInitialization, std::move(applyCutoffDuringInitialization));
}

inline void SettingsPopulator::addGaffAtomicChargesFile(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor atomicChargeFilePath("Path to atomic charges file for GAFF.");
  atomicChargeFilePath.setDefaultValue("");
  settings.push_back(SettingsNames::gaffAtomicChargesFile, std::move(atomicChargeFilePath));
}

inline void SettingsPopulator::addTitrationSiteFile(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor titrationSiteFile("Path to titation sites for MM parametrization.");
  titrationSiteFile.setDefaultValue("titrable_sites.dat");
  settings.push_back(SettingsNames::titrationSiteFile, std::move(titrationSiteFile));
}

inline void SettingsPopulator::addGaffAtomTypesFile(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor atomTypesFilePath("Path to atom types file for GAFF.");
  atomTypesFilePath.setDefaultValue("");
  settings.push_back(SettingsNames::gaffAtomTypesFile, std::move(atomTypesFilePath));
}

inline void SettingsPopulator::addQmAtomsOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntListDescriptor qmAtoms("A list containing the indices of the atoms in the QM region.");
  qmAtoms.setDefaultValue({});
  settings.push_back(Utils::SettingsNames::qmAtomsList, std::move(qmAtoms));
}

inline void SettingsPopulator::addChargeRedistributionOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor chargeRedistr(
      "Sets the charge redistribution scheme for the MM atoms close to the QM-MM boundary.");
  chargeRedistr.addOption(OptionNames::redistributedChargeOption);
  chargeRedistr.addOption(OptionNames::redistributedChargeAndDipolesOption);
  chargeRedistr.setDefaultOption(OptionNames::redistributedChargeOption);
  settings.push_back(SettingsNames::chargeRedistributionKey, std::move(chargeRedistr));
}

inline void SettingsPopulator::addElectrostaticEmbeddingOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor electrostaticEmbeddingOption(
      "Sets whether electrostatic embedding is used in QM/MM. The alternative is applying mechanical embedding only.");
  electrostaticEmbeddingOption.setDefaultValue(false);
  settings.push_back(Utils::SettingsNames::electrostaticEmbedding, std::move(electrostaticEmbeddingOption));
}

inline void SettingsPopulator::addQmRegionXyzFileOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor qmRegionXyzFile(
      "The path to a file to which the QM region can be dumped in XYZ format.");
  qmRegionXyzFile.setDefaultValue("");
  settings.push_back(SettingsNames::qmRegionXyzFile, std::move(qmRegionXyzFile));
}

inline void SettingsPopulator::addIgnoreQmOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor ignoreQm(
      "Whether to ignore all contributions from the QM calculation, and therefore, not performing it.");
  ignoreQm.setDefaultValue(false);
  settings.push_back(Utils::SettingsNames::ignoreQmOption, std::move(ignoreQm));
}

inline void SettingsPopulator::addOptimizeLinksOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor optLinks(
      "Whether to optimize the position of the link nuclei before reporting an energy.");
  optLinks.setDefaultValue(false);
  settings.push_back(Utils::SettingsNames::optimizeLinks, std::move(optLinks));
}

inline void SettingsPopulator::addSilenceOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor silence("Whether to silence the standard output of the subcalculators.");
  silence.setDefaultValue(false);
  settings.push_back(SettingsNames::silenceUnderlyingCalculators, std::move(silence));
}

inline void SettingsPopulator::addReducedQmMmEnergyOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor reducedEnergyOption(
      "Sets whether an additional MM calculation shall be performed to evaluate the reduced QM/MM energy without any "
      "MM contributions for atoms located solely within the environment.");
  reducedEnergyOption.setDefaultValue(false);
  settings.push_back(SettingsNames::calculateReducedQmMmEnergy, std::move(reducedEnergyOption));
}

inline void SettingsPopulator::addBondOrderThreshold(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor bondOrderThreshold(
      "Sets the threshold for which bond orders to consider as bonds.");
  bondOrderThreshold.setMinimum(0.0);
  bondOrderThreshold.setMaximum(2.0);
  bondOrderThreshold.setDefaultValue(0.5);
  settings.push_back(SettingsNames::bondOrderThreshold, std::move(bondOrderThreshold));
}

inline void SettingsPopulator::addConnectivityRefinementOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor refineConnectivity(
      "Whether the connectivities of the atoms obtained from covalent radii shall be refined applying quantum-chemical "
      "data.");
  refineConnectivity.setDefaultValue(true);
  settings.push_back(SettingsNames::refineConnectivity, std::move(refineConnectivity));
}

inline void SettingsPopulator::addExistingParameterFile(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor existingParameters(
      "Path to an MM parameter file, which already contains some parameters that can be re-used for the "
      "parametrization.");
  existingParameters.setDefaultValue("");
  settings.push_back(SettingsNames::existingParameters, std::move(existingParameters));
}

inline void SettingsPopulator::addConstrainMMParametersOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor constrainMMParameters(
      "Decides whether there should be constraints during the MM parameter optimization.");
  constrainMMParameters.setDefaultValue(true);
  settings.push_back(SettingsNames::constrainMMParameters, std::move(constrainMMParameters));
}

inline void
SettingsPopulator::addOptimizeImproperDihedralForceConstantsOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor optImproperDihedralForceConstants(
      "Decides whether the improper dihedral force constants of non-planar groups should be optimized during the MM "
      "parameter optimization.");
  optImproperDihedralForceConstants.setDefaultValue(true);
  settings.push_back(SettingsNames::optimizeImproperDihedralForceConstants, std::move(optImproperDihedralForceConstants));
}

inline void SettingsPopulator::addNumberAtomsThreshold(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor numberAtomsThreshold(
      "The maximum number of atoms where reference data is still generated for the whole system.");
  numberAtomsThreshold.setMinimum(1);
  numberAtomsThreshold.setMaximum(200);
  numberAtomsThreshold.setDefaultValue(150);
  settings.push_back(SettingsNames::numberAtomsThreshold, std::move(numberAtomsThreshold));
}

inline void SettingsPopulator::addGaussianMethodAndBasisSet(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor gaussianMethod(
      "The method used in the Gaussian calculation for the CM5 charges.");
  gaussianMethod.setDefaultValue("PBEPBE");
  settings.push_back(SettingsNames::gaussianMethod, std::move(gaussianMethod));

  Utils::UniversalSettings::StringDescriptor gaussianBasisSet(
      "The basis set used in the Gaussian calculation for the CM5 charges.");
  gaussianBasisSet.setDefaultValue("def2SVP");
  settings.push_back(SettingsNames::gaussianBasisSet, std::move(gaussianBasisSet));
}

inline void SettingsPopulator::addBaseWorkingDirectory(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor baseWorkingDirectory("Base directory for the calculations.");
  baseWorkingDirectory.setDefaultValue(Utils::FilesystemHelpers::currentDirectory());
  settings.push_back(Utils::ExternalQC::SettingsNames::baseWorkingDirectory, std::move(baseWorkingDirectory));
}

inline void SettingsPopulator::addExternalProgramNProcs(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor externalProgramNProcs(
      "Number of processes for a reference calculation by an external program.");
  externalProgramNProcs.setDefaultValue(1);
  externalProgramNProcs.setMinimum(1);
  settings.push_back(Utils::SettingsNames::externalProgramNProcs, std::move(externalProgramNProcs));
}

inline void SettingsPopulator::addReferenceMethodAndBasisSet(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor referenceMethod("The method used in reference calculations.");
  referenceMethod.setDefaultValue(""); // default will be set later depending on the program
  settings.push_back(SettingsNames::referenceMethod, std::move(referenceMethod));

  Utils::UniversalSettings::StringDescriptor referenceBasisSet("The basis set used in reference calculations.");
  referenceBasisSet.setDefaultValue(""); // default will be set later depending on the program
  settings.push_back(SettingsNames::referenceBasisSet, std::move(referenceBasisSet));
}

inline void SettingsPopulator::addSubsystemRadius(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor subsystemRadius(
      "The radius of the spheres defining the volume of the subsystems. Unit: Angstrom.");
  subsystemRadius.setMinimum(5.0);
  subsystemRadius.setMaximum(12.0);
  subsystemRadius.setDefaultValue(6.0);
  settings.push_back(SettingsNames::subsystemRadius, std::move(subsystemRadius));
}

inline void SettingsPopulator::addReferenceDataGenerationOptions(Utils::UniversalSettings::DescriptorCollection& settings,
                                                                 bool includeReadWriteMode) {
  if (includeReadWriteMode) {
    Utils::UniversalSettings::StringDescriptor refDataDir(
        "The name of the directory where the reference calculation files are loaded from and to which the molecular "
        "structure data is written.");
    refDataDir.setDefaultValue("reference_data");
    settings.push_back(SettingsNames::referenceDataDirectory, std::move(refDataDir));
  }

  Utils::UniversalSettings::OptionListDescriptor refDataMode(
      "Mode of the reference data generation (direct, reading, writing, database)");
  refDataMode.addOption(OptionNames::directMode);
  refDataMode.addOption(OptionNames::databaseMode);
  if (includeReadWriteMode) {
    refDataMode.addOption(OptionNames::readFromFilesMode);
    refDataMode.addOption(OptionNames::writeToFilesMode);
  }
  refDataMode.setDefaultOption(OptionNames::directMode);
  settings.push_back(SettingsNames::referenceDataMode, std::move(refDataMode));
}

inline void SettingsPopulator::addDatabaseSettings(Utils::UniversalSettings::DescriptorCollection& settings,
                                                   std::string defaultDatabaseName) {
  Utils::UniversalSettings::StringDescriptor databaseHost("The name or IP address of the database host.");
  databaseHost.setDefaultValue("localhost");
  settings.push_back(SettingsNames::databaseHost, std::move(databaseHost));

  Utils::UniversalSettings::StringDescriptor databaseName("The name of the database.");
  databaseName.setDefaultValue(defaultDatabaseName);
  settings.push_back(SettingsNames::databaseName, std::move(databaseName));

  Utils::UniversalSettings::IntDescriptor databasePort("The port through which to connect to the database.");
  databasePort.setMinimum(27000);
  databasePort.setMaximum(27999);
  databasePort.setDefaultValue(27017);
  settings.push_back(SettingsNames::databasePort, std::move(databasePort));
}

inline void SettingsPopulator::addDatabaseSleepTime(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor databaseSleepTime(
      "The sleep time in seconds inbetween database operations in reference data generation phase.");
  databaseSleepTime.setMinimum(1);
  databaseSleepTime.setMaximum(86400); // one day
  databaseSleepTime.setDefaultValue(60);
  settings.push_back(SettingsNames::databaseSleepTime, std::move(databaseSleepTime));
}

inline void SettingsPopulator::addAtomicInformationFile(Utils::UniversalSettings::DescriptorCollection& settings,
                                                        bool setEmptyDefault) {
  Utils::UniversalSettings::StringDescriptor atomicInfoFile(
      "Path to file containing the information about formal charges and unpaired electrons.");
  if (setEmptyDefault) {
    atomicInfoFile.setDefaultValue("");
  }
  else {
    atomicInfoFile.setDefaultValue("atomic_info.dat");
  }
  settings.push_back(SettingsNames::atomicInformationFile, std::move(atomicInfoFile));
}

inline void SettingsPopulator::addUseGaussianOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor useGaussianOption(
      "Decides whether to apply the Gaussian program to calculate CM5 charges.");
  useGaussianOption.setDefaultValue(false);
  settings.push_back(SettingsNames::useGaussianOptionKey, std::move(useGaussianOption));
}

inline void SettingsPopulator::addReferenceProgram(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::OptionListDescriptor referenceProgram(
      "Set program for the reference calculations in SFAM parametrization.");
  referenceProgram.addOption(OptionNames::turbomoleOption);
  referenceProgram.addOption(OptionNames::orcaOption);
  referenceProgram.addOption(OptionNames::sparrowOption);
  referenceProgram.addOption(OptionNames::xtbOption);
  referenceProgram.setDefaultOption(OptionNames::orcaOption);
  settings.push_back(SettingsNames::referenceProgram, std::move(referenceProgram));
}

inline void SettingsPopulator::addIncreaseScfSafetyOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor increaseScfSafety(
      "Decides whether options for a safer SCF convergence should be applied.");
  increaseScfSafety.setDefaultValue(false);
  settings.push_back(SettingsNames::increaseScfSafetyKey, std::move(increaseScfSafety));
}

inline void SettingsPopulator::addEarlyTerminationOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor earlyTermination(
      "Decides whether the reference data generation should be terminated early when enough data has been collected.");
  earlyTermination.setDefaultValue(true);
  settings.push_back(SettingsNames::enableEarlyTerminationKey, std::move(earlyTermination));
}

inline void SettingsPopulator::addReuseDatabaseOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor reuseDatabase("Decides whether an existing database should be employed as "
                                                         "the basis for the parametrization by reusing its data.");
  reuseDatabase.setDefaultValue(false);
  settings.push_back(SettingsNames::reuseDatabaseKey, std::move(reuseDatabase));
}

inline void
SettingsPopulator::addTerminateAfterReferenceDataGenerationOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor refDataOnly(
      "Decides whether to terminate program after reference data generation step of parametrization. Can only be set "
      "to true in case of the database or write mode.");
  refDataOnly.setDefaultValue(false);
  settings.push_back(SettingsNames::terminateAfterReferenceDataGeneration, std::move(refDataOnly));
}

inline void SettingsPopulator::addUseCsvInputFormatOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor useCsv("Whether to use the CSV input format during read mode.");
  useCsv.setDefaultValue(true);
  settings.push_back(SettingsNames::useCsvInputFormat, std::move(useCsv));
}

inline void SettingsPopulator::addConvertToCm5Option(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor convertChargesCm5(
      "Whether to convert atomic charges with the Charge Model 5 algorithm.");
  convertChargesCm5.setDefaultValue(true);
  settings.push_back(SettingsNames::convertChargesCm5, std::move(convertChargesCm5));
}

inline void SettingsPopulator::addYamlSettingsForDirectMode(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor yamlSettingsFilePath("");
  yamlSettingsFilePath.setDefaultValue("");
  settings.push_back(SettingsNames::yamlSettingsFilePath, std::move(yamlSettingsFilePath));
}

inline void SettingsPopulator::addQmRegionCenterAtoms(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntListDescriptor centerAtoms(
      "The center atom(s) around which the QM Region will be constructed.");
  centerAtoms.setDefaultValue({0});
  settings.push_back(SettingsNames::qmRegionCenterAtoms, std::move(centerAtoms));
}

inline void SettingsPopulator::addInitialRadius(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor initialRadius("The initial radius for the QM region generation.");
  initialRadius.setMinimum(2.0);
  initialRadius.setMaximum(12.0);
  initialRadius.setDefaultValue(6.0);
  settings.push_back(SettingsNames::initialRadiusForQmRegionSelection, std::move(initialRadius));
}

inline void SettingsPopulator::addCuttingProbability(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor cuttingProbability(
      "The cutting probability for the QM region generation.");
  cuttingProbability.setMinimum(0.1);
  cuttingProbability.setMaximum(1.0);
  cuttingProbability.setDefaultValue(1.0);
  settings.push_back(SettingsNames::cuttingProbability, std::move(cuttingProbability));
}

inline void SettingsPopulator::addQmRegionSizesSettings(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor minSize("Minimum size of the QM region.");
  minSize.setMinimum(1);
  minSize.setDefaultValue(100);
  settings.push_back(SettingsNames::qmRegionCandidateMinSize, std::move(minSize));

  Utils::UniversalSettings::IntDescriptor maxSize("Maximum size of the QM region.");
  maxSize.setMinimum(1);
  maxSize.setDefaultValue(120);
  settings.push_back(SettingsNames::qmRegionCandidateMaxSize, std::move(maxSize));

  Utils::UniversalSettings::IntDescriptor refMaxSize("Maximum size for the reference QM regions.");
  refMaxSize.setMinimum(1);
  refMaxSize.setDefaultValue(200);
  settings.push_back(SettingsNames::qmRegionRefMaxSize, std::move(refMaxSize));
}

inline void SettingsPopulator::addNumAttemptsPerRadiusOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor numAttemptsPerRadius(
      "Number of QM region generation attempts for a given initial radius.");
  numAttemptsPerRadius.setMinimum(1);
  numAttemptsPerRadius.setDefaultValue(100);
  settings.push_back(SettingsNames::numAttemptsPerRadius, std::move(numAttemptsPerRadius));
}

inline void SettingsPopulator::addMaxNumRefModelsOption(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor maxNumRefModels("Maximum number of reference models to be evaluated.");
  maxNumRefModels.setMinimum(1);
  maxNumRefModels.setDefaultValue(10);
  settings.push_back(SettingsNames::maxNumRefModels, std::move(maxNumRefModels));
}

inline void SettingsPopulator::addTolerancesForQmRegionSelection(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor tolPercentageError(
      "The tolerance in percent for the mean absolute error of forces for the QM region selection.");
  tolPercentageError.setMinimum(0.0);
  tolPercentageError.setDefaultValue(20.0);
  settings.push_back(SettingsNames::tolerancePercentageError, std::move(tolPercentageError));

  Utils::UniversalSettings::DoubleDescriptor tolPercentageSymmetryScore(
      "The tolerance in percent for the symmetry score based on the minimum symmetry score of all candidate models.");
  tolPercentageSymmetryScore.setMinimum(0.0);
  tolPercentageError.setDefaultValue(50.0);
  settings.push_back(SettingsNames::tolerancePercentageSymmetryScore, std::move(tolPercentageSymmetryScore));
}

inline void SettingsPopulator::addQmRegionSelectionRandomSeed(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor randomSeed("Random seed for the QM region selection.");
  randomSeed.setDefaultValue(42);
  settings.push_back(SettingsNames::qmRegionSelectionRandomSeed, std::move(randomSeed));
}

inline void SettingsPopulator::addPreparationDataDirectory(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor preparationDataDirectory("Base directory for the preparation.");
  preparationDataDirectory.setDefaultValue("preparation_data");
  settings.push_back(SettingsNames::preparationDataDirectory, std::move(preparationDataDirectory));
}

inline void SettingsPopulator::addSolvation(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor solvateStructure("Solvates the input structure.");
  solvateStructure.setDefaultValue(false);
  settings.push_back(SettingsNames::solvateStructure, std::move(solvateStructure));
}

inline void SettingsPopulator::addNumberOfSolventShells(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::IntDescriptor numberOfSolventShells("Defines the number of solvent shells");
  numberOfSolventShells.setDefaultValue(1);
  settings.push_back(SettingsNames::numberOfSolventShells, std::move(numberOfSolventShells));
}

inline void SettingsPopulator::addChargedCandNTermini(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor chargedCandNTermini(
      "Whether to protonate C and N termini in their charged (NH3+, COO-) form");
  chargedCandNTermini.setDefaultValue(true);
  settings.push_back(SettingsNames::chargedCandNTermini, std::move(chargedCandNTermini));
}

inline void SettingsPopulator::addPhValueOfSystem(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::DoubleDescriptor phValueOfSystem("pH value for the protonation.");
  phValueOfSystem.setDefaultValue(7.0);
  settings.push_back(SettingsNames::phValueOfSystem, std::move(phValueOfSystem));
}

inline void SettingsPopulator::addTitration(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor titrate(
      "Titrates the structure to protonate amino acids according to the pH value.");
  titrate.setDefaultValue(false);
  settings.push_back(SettingsNames::titrate, std::move(titrate));
}

inline void SettingsPopulator::addUseThermochemistryForTitration(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::BoolDescriptor useThermoChemistryForTitration(
      "Sets whether thermochemical data (gibbs free energy) should be used for calculation of pka values and "
      "protonation states.");
  useThermoChemistryForTitration.setDefaultValue(false);
  settings.push_back(SettingsNames::useThermoChemistryForTitration, std::move(useThermoChemistryForTitration));
}

inline void SettingsPopulator::addTrainingDataDirectory(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor trainingDataDirectory(
      "Sets the directory in which the reference calculation for the training data are stored.");
  trainingDataDirectory.setDefaultValue("");
  settings.push_back(SettingsNames::trainingDataDirectory, std::move(trainingDataDirectory));
}

inline void SettingsPopulator::addMethodFamily(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor methodFamily("The method family such as 'DFT/SFAM'.");
  methodFamily.setDefaultValue("PM6/SFAM");
  settings.push_back(Utils::SettingsNames::methodFamily, std::move(methodFamily));
}

inline void SettingsPopulator::addProgram(Utils::UniversalSettings::DescriptorCollection& settings) {
  Utils::UniversalSettings::StringDescriptor program("The underlying programs such as 'Turbomole/SFAM'.");
  program.setDefaultValue("Any/Swoose");
  settings.push_back(Utils::SettingsNames::program, std::move(program));
}

} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_SETTINGSPOPULATOR_H
