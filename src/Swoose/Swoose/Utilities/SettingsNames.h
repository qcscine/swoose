/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_SETTINGSNAMES_H
#define SWOOSEUTILITIES_SETTINGSNAMES_H

#include <Utils/Settings.h>

namespace Scine {
namespace SwooseUtilities {
namespace SettingsNames {

// Mostly used for MM
static constexpr const char* onlyCalculateBondedContribution = "covalent_contributions_only";
static constexpr const char* printContributionsMolecularMechanics = "print_mm_contributions";
static constexpr const char* sfamAtomTypeLevel = "atom_type_level";
static constexpr const char* connectivityFilePath = "mm_connectivity_file";
static constexpr const char* detectBondsWithCovalentRadii = "covalent_radii_bond_detection";
static constexpr const char* nonCovalentCutoffRadius = "non_covalent_cutoff";
static constexpr const char* hydrogenBondCorrection = "hydrogen_bond_correction";
static constexpr const char* applyCutoffDuringInitialization = "apply_cutoff_during_initialization";
static constexpr const char* gaffAtomicChargesFile = "gaff_atomic_charges_file";
static constexpr const char* gaffAtomTypesFile = "gaff_atom_types_file";

// Mostly used for QM/MM
static constexpr const char* chargeRedistributionKey = "charge_redistribution";
static constexpr const char* calculateReducedQmMmEnergy = "reduced_qmmm_energy";
static constexpr const char* qmRegionXyzFile = "qm_region_file";
static constexpr const char* silenceUnderlyingCalculators = "silence_underlying_calculators";

// Mostly used for MM parametrization
static constexpr const char* bondOrderThreshold = "bond_order_threshold";
static constexpr const char* refineConnectivity = "refine_connectivity_qm";
static constexpr const char* existingParameters = "existing_parameters";
static constexpr const char* constrainMMParameters = "constrain_mm_parameters";
static constexpr const char* optimizeImproperDihedralForceConstants = "optimize_improper_dihedral_force_constants";
static constexpr const char* numberAtomsThreshold = "number_atoms_threshold";
static constexpr const char* subsystemRadius = "subsystem_radius";
static constexpr const char* referenceDataDirectory = "ref_data_directory";
static constexpr const char* referenceDataMode = "ref_data_mode";
static constexpr const char* databaseHost = "database_host";
static constexpr const char* databaseName = "database_name";
static constexpr const char* databasePort = "database_port";
static constexpr const char* databaseSleepTime = "database_sleep_time";
static constexpr const char* atomicInformationFile = "atomic_info_file";
static constexpr const char* referenceMethod = "reference_method";
static constexpr const char* referenceBasisSet = "reference_basis_set";
static constexpr const char* referenceProgram = "reference_program";
static constexpr const char* gaussianMethod = "gaussian_method";
static constexpr const char* gaussianBasisSet = "gaussian_basis_set";
static constexpr const char* useGaussianOptionKey = "use_gaussian";
static constexpr const char* increaseScfSafetyKey = "increase_scf_safety";
static constexpr const char* enableEarlyTerminationKey = "enable_early_termination";
static constexpr const char* reuseDatabaseKey = "reuse_database";
static constexpr const char* terminateAfterReferenceDataGeneration = "ref_data_generation_only";
static constexpr const char* useCsvInputFormat = "use_csv";
static constexpr const char* convertChargesCm5 = "convert_charges_cm5";
static constexpr const char* yamlSettingsFilePath = "yaml_settings_file_path";
static constexpr const char* titrate = "titrate";
static constexpr const char* useThermoChemistryForTitration = "use_thermochemistry_for_titration";
static constexpr const char* trainingDataDirectory = "training_data_directory";
static constexpr const char* titrationSiteFile = "titration_site_file";

// Mostly used for QM/MM region selection
static constexpr const char* qmRegionCenterAtoms = "qm_region_center_atoms";
static constexpr const char* initialRadiusForQmRegionSelection = "initial_radius";
static constexpr const char* cuttingProbability = "cutting_probability";
static constexpr const char* qmRegionCandidateMinSize = "qm_region_min_size";
static constexpr const char* qmRegionCandidateMaxSize = "qm_region_max_size";
static constexpr const char* qmRegionRefMaxSize = "ref_max_size";
static constexpr const char* numAttemptsPerRadius = "num_attempts_per_radius";
static constexpr const char* maxNumRefModels = "max_num_ref_models";
static constexpr const char* tolerancePercentageError = "tol_percentage_error";
static constexpr const char* tolerancePercentageSymmetryScore = "tol_percentage_sym_score";
static constexpr const char* qmRegionSelectionRandomSeed = "qm_region_selection_seed";

// Mostly used for structure preparation
static constexpr const char* preparationDataDirectory = "preparation_directory";
static constexpr const char* atomicInfoFile = "atomic_info_file";
static constexpr const char* solvateStructure = "solvate_structure";
static constexpr const char* numberOfSolventShells = "num_solvent_shells";
static constexpr const char* phValueOfSystem = "pH_value_for_protonation";
static constexpr const char* chargedCandNTermini = "charged_termini";

} // namespace SettingsNames
} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_SETTINGSNAMES_H
