/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Parametrizer.h"
#include "CalculationManager.h"
#include "MMParametrizationSettings.h"
#include "MolecularSystemPartitioner.h"
#include "OptimizationSetup.h"
#include "ParameterOptimizer.h"
#include "ParametrizationUtils/AtomicChargesAssembler.h"
#include "ParametrizationUtils/ConnectivityGenerator.h"
#include "ParametrizationUtils/FullHessianAssembler.h"
#include "ParametrizationUtils/ParameterFileWriter.h"
#include "ParametrizationUtils/ReparametrizationHelper.h"
#include "ParametrizationUtils/SuperfluousFragmentIdentifier.h"
#include <Core/Log.h>
#include <Swoose/MolecularMechanics/SFAM/SfamAtomTypeIdentifier.h>
#include <Swoose/MolecularMechanics/Topology/IndexedStructuralTopologyCreator.h>
#include <Swoose/Utilities/AtomicInformationReader.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Utils/Geometry/ElementInfo.h>
#include <chrono>

namespace Scine {
namespace MMParametrization {

Parametrizer::Parametrizer() {
  this->settings_ = std::make_shared<MMParametrizationSettings>();
  this->connectivityGenerator_ = std::make_unique<ConnectivityGenerator>(this->data_, this->settings_, this->getLog());
}

std::string Parametrizer::name() const {
  return "SFAM_parametrizer";
}

// Performs the parametrization in 5 steps
void Parametrizer::parametrize(Utils::AtomCollection structure) {
  auto start =
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  performInitialSetup(std::move(structure));

  auto t1 =
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  generateReferenceData();
  // If the data generation mode is "write" it does not make any sense to continue from here
  if (settings_->getString(SwooseUtilities::SettingsNames::referenceDataMode) == SwooseUtilities::OptionNames::writeToFilesMode)
    return;
  // Stop at this point if it was requested
  if (settings_->getBool(SwooseUtilities::SettingsNames::terminateAfterReferenceDataGeneration)) {
    this->getLog().output << Core::Log::nl
                          << "The program terminates now, because the reference data generation is done. This was "
                             "requested via the settings."
                          << Core::Log::endl;
    return;
  }

  auto t2 =
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  setupParameterOptimization();

  auto t3 =
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  optimizeParameters();

  auto t4 =
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  writeParametersAndConnectivity();

  auto end =
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  // Print out timings
  this->getLog().output << "Parametrization procedure completed. Timings:" << Core::Log::nl << Core::Log::nl;
  this->getLog().output << "Initial set up and fragmentation: " << (t1 - start) / 1000.0 << " seconds." << Core::Log::nl;
  this->getLog().output << "Reference data generation: " << (t2 - t1) / 1000.0 << " seconds." << Core::Log::nl;
  this->getLog().output << "Set up of optimization: " << (t3 - t2) / 1000.0 << " seconds." << Core::Log::nl;
  this->getLog().output << "Parameter optimization: " << (t4 - t3) / 1000.0 << " seconds." << Core::Log::nl;
  this->getLog().output << "Writing files: " << (end - t4) / 1000.0 << " seconds." << Core::Log::nl;
  this->getLog().output << "Total time: " << (end - start) / 1000.0 << " seconds." << Core::Log::endl;
}

void Parametrizer::performInitialSetup(Utils::AtomCollection structure) {
  settings_->normalizeStringCases(); // convert all option names to lower case letters
  if (!settings_->valid()) {
    settings_->throwIncorrectSettings();
  }

  // Set defaults for reference_method and reference_basis_set settings if not set by user
  bool referenceMethodWasNotSelected = settings_->getString(SwooseUtilities::SettingsNames::referenceMethod).empty();
  if (referenceMethodWasNotSelected)
    setDefaultsForMethodAndBasisSetSettings();

  // Move structure to data object
  data_.fullStructure = std::move(structure);
  // Set number of atoms in data object
  data_.numberOfAtoms = data_.fullStructure.size();

  performAdditionalSettingsChecks();

  this->getLog().output << "Starting the parametrization of the MM model..." << Core::Log::endl;

  // Read the atomic information file containing information about the formal charges and unpaired electrons
  std::string atomicInfoFile = settings_->getString(SwooseUtilities::SettingsNames::atomicInformationFile);
  if (!atomicInfoFile.empty()) {
    SwooseUtilities::AtomicInformationReader atomicInfoReader(this->getLog());
    atomicInfoReader.read(atomicInfoFile, data_.formalCharges, data_.unpairedElectrons, data_.numberOfAtoms);
  }

  // Generate connectivity and update the listsOfNeighbors member of data_
  connectivityGenerator_->generateInitialListsOfNeighbors();
  // Generate topology and atom types
  generateTopology();
  generateAtomTypes();

  // If a parameter file of already existing parameters is given, create helper class for re-parametrization strategy
  auto existingParametersPath = settings_->getString(SwooseUtilities::SettingsNames::existingParameters);
  if (!existingParametersPath.empty()) {
    reparametrizationHelper_ = std::make_shared<ReparametrizationHelper>(data_, this->getLog());
    reparametrizationHelper_->parseProvidedParameters(existingParametersPath);
    reparametrizationHelper_->manipulateTopology();
  }

  if (data_.numberOfAtoms > settings_->getInt(SwooseUtilities::SettingsNames::numberAtomsThreshold))
    this->getLog().output << "Now dividing the system into its fragments..." << Core::Log::endl;
  // Partition system into subsystems
  MolecularSystemPartitioner partitioner(data_, settings_, this->getLog());
  partitioner.divideSystem();
  // Identify superfluous fragments
  SuperfluousFragmentIdentifier::identifySuperfluousFragments(data_, this->getLog(), reparametrizationHelper_);
}

void Parametrizer::generateReferenceData() {
  CalculationManager calculationManager(data_, settings_, this->getLog());
  calculationManager.calculateReferenceData();
}

void Parametrizer::setupParameterOptimization() {
  // Refine the connectivity based on quantum-chemically calculated bond orders if required
  if (settings_->getBool(SwooseUtilities::SettingsNames::refineConnectivity)) {
    connectivityGenerator_->refineListsOfNeighbors();
    // Generate topology and atom types based on refined connectivity
    generateTopology();
    generateAtomTypes();
    if (reparametrizationHelper_ != nullptr)
      reparametrizationHelper_->manipulateTopology();
  }

  // Assemble the atomic charges for the whole system in the data object from the subsystem results.
  AtomicChargesAssembler::assembleAtomicCharges(data_);
  AtomicChargesAssembler::renormalizeAtomicCharges(data_);
  std::vector<std::vector<double>>().swap(data_.atomicChargesForEachFragment); // Free memory

  // Assemble full Hessian matrix from subsystem Hessians
  FullHessianAssembler hessianAssembler(data_, this->getLog());
  hessianAssembler.assembleFullHessian();

  // Free memory
  data_.vectorOfStructures.clear();

  // Generate initial parameters
  OptimizationSetup optSetup(data_, settings_);
  optSetup.generateInitialParameters();

  // Free more memory
  data_.vectorOfHessians.clear();
  data_.vectorOfOptimizedStructures.clear();
}

void Parametrizer::optimizeParameters() {
  ParameterOptimizer optimizer(data_, settings_, this->getLog());
  optimizer.optimizeParameters();
}

void Parametrizer::writeParametersAndConnectivity() {
  std::string parameterFilePath = settings_->getString(SwooseUtilities::SettingsNames::parameterFilePath);
  std::string connectivityFilePath = settings_->getString(SwooseUtilities::SettingsNames::connectivityFilePath);
  if (!parameterFilePath.empty())
    ParameterFileWriter::writeSfamParametersToFile(parameterFilePath, data_.parameters);
  if (!connectivityFilePath.empty())
    SwooseUtilities::ConnectivityFileHandler::writeListsOfNeighbors(connectivityFilePath, data_.listsOfNeighbors);
}

void Parametrizer::generateTopology() {
  MolecularMechanics::IndexedStructuralTopologyCreator topologyCreator(data_.listsOfNeighbors);
  data_.topology = topologyCreator.calculateIndexedStructuralTopology();

  topologyCreator.addHydrogenBondsToIndexedStructuralTopology(data_.topology, data_.fullStructure);
}

void Parametrizer::generateAtomTypes() {
  MolecularMechanics::SfamAtomTypeIdentifier atomTypeIdentifier(data_.numberOfAtoms, data_.fullStructure.getElements(),
                                                                data_.listsOfNeighbors);
  MolecularMechanics::SfamAtomTypeLevel atl = MolecularMechanics::SfamAtomTypeIdentifier::generateSfamAtomTypeLevelFromString(
      settings_->getString(SwooseUtilities::SettingsNames::sfamAtomTypeLevel));
  data_.atomTypes = atomTypeIdentifier.getAtomTypes(atl);
}

void Parametrizer::performAdditionalSettingsChecks() {
  auto refDataMode = settings_->getString(SwooseUtilities::SettingsNames::referenceDataMode);

#ifndef SWOOSE_COMPILE_DATABASE
  if (refDataMode == SwooseUtilities::OptionNames::databaseMode)
    throw std::runtime_error("Swoose was not compiled with database support.");
#endif

  // Check whether the direct mode for reference calculations is valid.
  if (refDataMode == SwooseUtilities::OptionNames::directMode) {
    if (data_.numberOfAtoms > settings_->getInt(SwooseUtilities::SettingsNames::numberAtomsThreshold)) {
      throw std::runtime_error(
          "The direct mode for reference calculations cannot be applied if the molecular system has to be fragmented.");
    }
  }
  // Check whether the 'ref_data_generation_only' option is allowed.
  if (settings_->getBool(SwooseUtilities::SettingsNames::terminateAfterReferenceDataGeneration)) {
    if (refDataMode != SwooseUtilities::OptionNames::databaseMode &&
        refDataMode != SwooseUtilities::OptionNames::writeToFilesMode) {
      throw std::runtime_error("The setting \"" +
                               std::string(SwooseUtilities::SettingsNames::terminateAfterReferenceDataGeneration) +
                               "\" cannot be set in read and direct mode.");
    }
  }

  auto referenceProgram = settings_->getString(SwooseUtilities::SettingsNames::referenceProgram);
  if (referenceProgram == SwooseUtilities::OptionNames::xtbOption ||
      referenceProgram == SwooseUtilities::OptionNames::sparrowOption) {
    if (settings_->getBool(SwooseUtilities::SettingsNames::convertChargesCm5))
      this->getLog().warning << "Warning: With xTB or Sparrow as the reference program, we recommend to set '"
                             << SwooseUtilities::SettingsNames::convertChargesCm5 << "' to false." << Core::Log::nl
                             << Core::Log::endl;
  }
}

void Parametrizer::setDefaultsForMethodAndBasisSetSettings() {
  auto selectedProgram = settings_->getString(SwooseUtilities::SettingsNames::referenceProgram);
  std::string defaultReferenceMethodString = "PBE D3BJ";   // ORCA option base case
  std::string defaultReferenceBasisSetString = "def2-SVP"; // ORCA option base case
  if (selectedProgram == SwooseUtilities::OptionNames::xtbOption) {
    defaultReferenceMethodString = "gfn2";
    defaultReferenceBasisSetString = "";
  }
  else if (selectedProgram == SwooseUtilities::OptionNames::sparrowOption) {
    defaultReferenceMethodString = "pm6";
    defaultReferenceBasisSetString = "";
  }
  else if (selectedProgram == SwooseUtilities::OptionNames::turbomoleOption) {
    defaultReferenceMethodString = "pbe D3BJ"; // lower-case functional name
    defaultReferenceBasisSetString = "def2-SVP";
  }
  settings_->modifyString(SwooseUtilities::SettingsNames::referenceMethod, defaultReferenceMethodString);
  settings_->modifyString(SwooseUtilities::SettingsNames::referenceBasisSet, defaultReferenceBasisSetString);
}

const Utils::Settings& Parametrizer::settings() const {
  return *settings_;
}

Utils::Settings& Parametrizer::settings() {
  return *settings_;
}

} // namespace MMParametrization
} // namespace Scine
