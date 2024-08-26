/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "StructureProcessor.h"
#include "ProteinStructures.h"
#include "Protonation/ProtonationHandler.h"
#include "SpecialCaseHandler.h"
#include "StructurePreparationHelper.h"
#include "StructurePreparationIO.h"
#include "StructurePreparationSettings.h"
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/TitrationFileHandler.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/ChemicalFileFormats/PdbStreamHandler.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/IO/NativeFilenames.h>
#include <numeric>

namespace Scine {
namespace StructurePreparation {

StructureProcessor::StructureProcessor() {
  this->settings_ = std::make_shared<StructurePreparationSettings>();
}

// Performs the PDB Preparation in 3 steps
void StructureProcessor::prepare(const std::string& structureFile, const std::string& mode) {
  performInitialSetup(structureFile);

  if (mode == SwooseUtilities::OptionNames::analyzeMode) {
    analyzeStructure(structureFile);

    log_.output << "Done. " << Core::Log::nl;
    if (boost::filesystem::exists(files_.nonRegContainerFile)) {
      log_.output << "Please look at " + files_.nonRegContainerFile + " and make changes accordingly. "
                  << "\n"
                  << "You can then restart the pdb preparation app with preparation_mode: "
                  << SwooseUtilities::OptionNames::protonationMode << Core::Log::nl;
    }
  }
  else if (mode == SwooseUtilities::OptionNames::protonationMode) {
    protonate(); // protonate entire structure in uncharged reference state

    if (boost::filesystem::exists(files_.protonatedNonRegContainerFile)) {
      log_.output << "Please check the protonation at " + files_.protonatedNonRegContainerFile << "\n"
                  << "You can then restart the PDB preparation app with preparation_mode "
                  << SwooseUtilities::OptionNames::finalizeMode << Core::Log::nl;
    }
  }
  else if (mode == SwooseUtilities::OptionNames::finalizeMode) {
    finalize();
  }
  else if (mode == SwooseUtilities::OptionNames::automateMode) {
    StructurePreparationHelper::performInitialCheck(files_, 1, settings_);
    analyzeStructure(structureFile);
    if (!protonationDone_) {
      protonate();
    }
    if (!finalizationDone_) {
      finalize();
    }
  }
}

void StructureProcessor::analyzeStructure(const std::string& structureFile) {
  log_.output << "Analyzing the input structure ..." << Core::Log::nl;
  StructurePreparationHelper::performInitialCheck(files_, 1, settings_);
  auto structure = getAtomCollectionFromInput(structureFile, false);
  structure.removeAtomsByResidueLabel({"HOH"});
  StructurePreparationHelper::updatePdbPreparationData(data_, structure);
  StructurePreparationHelper::performGraphAnalysisOnStructure(data_, data_.fullStructure);
  StructurePreparationHelper::updateNonRegContainerVector(data_);
  StructurePreparationHelper::reevaluateConnectivityForAminoAcids(data_);
  if (size_t(data_.vectorOfNonRegContainerIndices.size() + data_.vectorOfProteinIndices.size()) !=
      size_t(data_.fullStructure.size())) {
    log_.output << "Structure analysis failed!" << Core::Log::nl;
  }
  StructurePreparationIO::writePdbFileWithResidueSpecifier(data_, files_.proteinFile, log_);
  processNonRegContainerSubstructure();
}

void StructureProcessor::protonate() {
  log_.output << "Protonating protein (and nonRegContainer, if present) ... " << Core::Log::nl;
  StructurePreparationHelper::performInitialCheck(files_, 2, settings_);
  ProtonationHandler protonationHandler(data_, files_, settings_);
  if (boost::filesystem::exists(files_.proteinFile)) {
    auto structure = Utils::ChemicalFileHandler::read(files_.proteinFile).first;
    Utils::ChemicalFileHandler::write(
        Utils::NativeFilenames::combinePathSegments(files_.workingDirectory, "protein.xyz"), structure);
    protonationHandler.setProteinStructure(structure);
    protonationHandler.protonateAllAminoAcids();
    // writeAtomicInfoFileForProtein();
  }
  if (boost::filesystem::exists(files_.nonRegContainerFile)) {
    protonationHandler.protonateNonRegContainerWithExternalSoftware();
  }

  log_.output << "Protonation done." << Core::Log::nl;
}

void StructureProcessor::finalize() {
  log_.output << "Finalizing ... " << Core::Log::nl;
  StructurePreparationHelper::performInitialCheck(files_, 3, settings_);
  // merge substructures
  StructurePreparationHelper::mergeProteinAndNonRegContainer(data_, files_);
  // handle the boundaries
  StructurePreparationHelper::handleBoundariesBetweenProteinAndNonRegContainer(data_);
  // write the atomic info file for the protein first
  StructurePreparationIO::writeAtomicInfoFileForProtein(data_, files_.atomicInfoFile);
  StructurePreparationIO::addAtomicInformationForNonRegContainer(files_, data_.subsystemMapping);
  auto titrableSites = StructurePreparationHelper::collectTitrableSites(data_);
  if (!titrableSites.empty()) {
    SwooseUtilities::TitrationFileHandler::writeTitrationSitesFile(files_.titrationSitesFile, titrableSites);
    for (const auto& site : titrableSites) {
      log_.debug << site.residueName << " : " << site.criticalAtom << Core::Log::nl;
    }
  }

  if (settings_->getBool(SwooseUtilities::SettingsNames::solvateStructure)) {
    log_.output << "Solvating the structure ..." << Core::Log::nl;

    auto solvent = StructurePreparationHelper::addSolvation(data_, settings_);

    Utils::ChemicalFileHandler::write(
        Utils::NativeFilenames::combinePathSegments(files_.workingDirectory, "solvent.xyz"), solvent);
    // Merge solvent and solute
    Utils::AtomCollection totalSystem;
    totalSystem += solvent;
    totalSystem += data_.fullStructure;
    StructurePreparationHelper::updatePdbPreparationData(data_, totalSystem);
  }
  SwooseUtilities::ConnectivityFileHandler::writeListsOfNeighbors(files_.connectivityFile, data_.listsOfNeighbors);
  Utils::ChemicalFileHandler::write(files_.systemFile, data_.fullStructure);
  log_.output << "Finalization done. You can see the final structure at " + files_.systemFile << Core::Log::nl;
}

void StructureProcessor::performInitialSetup(const std::string& structureFile) {
  files_.preparationDataDirectory = settings_->getString(SwooseUtilities::SettingsNames::preparationDataDirectory);

  // Create working directory
  boost::filesystem::path path(structureFile);
  std::string systemName =
      Utils::NativeFilenames::removeExtension(Utils::NativeFilenames::removeTrailingSeparator(path.filename().string()));
  files_.workingDirectory = Utils::NativeFilenames::combinePathSegments(files_.preparationDataDirectory, systemName);
  Utils::FilesystemHelpers::createDirectories(files_.workingDirectory);

  files_.initialize();

  files_.atomicInfoFile = Utils::NativeFilenames::combinePathSegments(
      files_.workingDirectory, settings_->getString(SwooseUtilities::SettingsNames::atomicInfoFile));
  // write connectivity file
  files_.connectivityFile = Utils::NativeFilenames::combinePathSegments(
      files_.workingDirectory, settings_->getString(SwooseUtilities::SettingsNames::connectivityFilePath));
}

Utils::AtomCollection StructureProcessor::getAtomCollectionFromInput(const std::string& structureFile, bool includeH) const {
  std::string suffix = StructurePreparationIO::getSuffix(structureFile);
  bool isPdb = suffix == "pdb";
  // Save structure as AtomCollection
  if (isPdb) {
    Utils::PdbStreamHandler handler;
    handler.setReadH(includeH);
    // first read in all solvent molecules
    std::ifstream input;
    input.open(structureFile);
    auto structures = handler.read(input);

    input.clear();
    input.seekg(0);

    std::vector<Utils::AtomCollection> solvents;
    try {
      solvents = handler.read(input);
      for (auto& solvent : solvents) {
        solvent.keepAtomsByResidueLabel({"HOH"});
      }
    }
    catch (...) {
      solvents.resize(0);
    }

    if (structures.size() > 1) {
      int counter = 0;
      // Write all overlaying structures to separate files
      for (const auto& structure : structures) {
        std::string structureFileName = "structure_" + std::to_string(counter) + ".xyz";
        std::string solventFileName = "solvent_" + std::to_string(counter) + ".xyz";
        Utils::ChemicalFileHandler::write(
            Utils::NativeFilenames::combinePathSegments(files_.workingDirectory, structureFileName), structure);
        if (solvents.size() > 1 && (solvents.size() == structures.size())) {
          Utils::ChemicalFileHandler::write(
              Utils::NativeFilenames::combinePathSegments(files_.workingDirectory, solventFileName), solvents.at(counter));
        }
        else if (solvents.size() == 1) {
          Utils::ChemicalFileHandler::write(
              Utils::NativeFilenames::combinePathSegments(files_.workingDirectory, solventFileName), solvents.at(0));
        }
        counter++;
      }
    }
    if (structures[0].size() == 0)
      throw std::runtime_error("Structure could not be parsed from the PDB file. ");
    return structures[0];
  }

  else {
    Utils::AtomCollection xyzStructure = Utils::ChemicalFileHandler::read(structureFile).first;
    if (includeH) {
      return xyzStructure;
    }
    Utils::AtomCollection structure;
    for (int i = 0; i < xyzStructure.size(); ++i) {
      // omit hydrogens
      if (xyzStructure.getElement(i) != Utils::ElementType::H)
        structure.push_back(xyzStructure.at(i));
    }
    return structure;
  }
}

void StructureProcessor::processNonRegContainerSubstructure() {
  int numberOfNonRegContainerAtoms = data_.numberOfAtoms - data_.vectorOfProteinIndices.size();

  if (numberOfNonRegContainerAtoms == 0) {
    log_.output << "No non-amino-acid residues were detected in your input structure. Directly proceeding to "
                   "protonation and finalization .. "
                << Core::Log::nl;
    protonate();
    protonationDone_ = true;
    finalize();
    finalizationDone_ = true;
  }
  else {
    Utils::AtomCollection nonRegContainer(numberOfNonRegContainerAtoms);

    int counter = 0;
    for (auto index : data_.vectorOfNonRegContainerIndices) {
      Utils::ElementType element = data_.fullStructure.getElement(index);
      Utils::Position position = data_.fullStructure.getPosition(index);
      nonRegContainer.setElement(counter, element);
      nonRegContainer.setPosition(counter, position);
      counter++;
    }
    Utils::ChemicalFileHandler::write(files_.nonRegContainerFile, nonRegContainer);
  }
}

const Utils::Settings& StructureProcessor::settings() const {
  return *settings_;
}

Utils::Settings& StructureProcessor::settings() {
  return *settings_;
}

void StructureProcessor::setLog(Core::Log& log) {
  log_ = std::move(log);
}

std::string StructureProcessor::name() const {
  return "PDB_Preparator";
}

void StructureProcessor::setFiles(StructurePreparationFiles files) {
  files_ = std::move(files);
}

} // namespace StructurePreparation
} // namespace Scine
