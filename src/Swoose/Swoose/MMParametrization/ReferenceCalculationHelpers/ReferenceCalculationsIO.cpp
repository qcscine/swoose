/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ReferenceCalculationsIO.h"
#include "../ParametrizationData.h"
#include <Core/Log.h>
#include <Swoose/MMParametrization/MMParametrizationSettings.h>
#include <Swoose/MMParametrization/ParametrizationUtils/TitrationHelper.h>
#include <Swoose/StructurePreparation/ProteinStructures.h>
#include <Swoose/StructurePreparation/Protonation/TitrationData.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/FragmentationHelper.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/ExternalQC/Gaussian/GaussianOutputParser.h>
#include <Utils/ExternalQC/Orca/OrcaHessianOutputParser.h>
#include <Utils/ExternalQC/Orca/OrcaMainOutputParser.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleMainOutputParser.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Properties/AtomicCharges/ChargeModel5.h>
#include <boost/filesystem.hpp>
#include <fstream>

namespace Scine {
namespace MMParametrization {
namespace ReferenceCalculationsIO {

void writeXyzFiles(ParametrizationData& data, std::string referenceDataDir, std::shared_ptr<Utils::Settings> settings) {
  // Write the structures to disk as xyz files along with information about charge, multiplicity and constrained atoms.
  for (int i = 0; i < int(data.vectorOfStructures.size()); ++i) {
    Utils::FilesystemHelpers::createDirectories(
        Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i)));

    // Write nothing if it is a superfluous fragment.
    if (std::find(data.superfluousFragments.begin(), data.superfluousFragments.end(), i) != data.superfluousFragments.end())
      continue;

    // Write structure
    std::string filename;
    if (data.vectorOfStructures.size() == 1)
      filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "molecule.xyz");
    else
      filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "fragment.xyz");
    Utils::ChemicalFileHandler::write(filename, *data.vectorOfStructures[i]);

    // Write charge and spin multiplicity
    std::string infoFilename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "info.dat");
    std::ofstream infoFile(infoFilename);
    infoFile << data.vectorOfChargesAndMultiplicities[i].first << "  "
             << data.vectorOfChargesAndMultiplicities[i].second << std::endl;
    infoFile.close();

    // Write constrained atoms
    auto constrainedAtoms = data.constrainedAtoms[i];
    std::string constrFilename =
        Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "constr.dat");
    writeConstrainedAtomsFile(constrainedAtoms, constrFilename);

    // Write additional data for titration
    if (settings->getBool(SwooseUtilities::SettingsNames::titrate)) {
      if (data.vectorOfStructures.size() > 1 && data.pHSensitiveSites.count(i) == 1)
        writeAdditionalDataForTitration(data, i, 0, referenceDataDir, settings);
      else if (data.vectorOfStructures.size() == 1 && data.pHSensitiveSites.size() == 1)
        writeAdditionalDataForTitration(data, i, data.pHSensitiveSites.begin()->first, referenceDataDir, settings);
    }
  }
}

void writeAdditionalDataForTitration(ParametrizationData& data, int fragmentIndex, int criticalAtomIndex,
                                     std::string referenceDataDir, std::shared_ptr<Utils::Settings> settings) {
  std::string aminoAcid;
  if (data.vectorOfStructures.size() == 1)
    aminoAcid = data.pHSensitiveSites.find(criticalAtomIndex)->second;
  else
    aminoAcid = data.pHSensitiveSites.find(fragmentIndex)->second;
  StructurePreparation::AminoAcidCategorizer categorizer;
  bool isBase =
      (std::find(categorizer.bases.begin(), categorizer.bases.end(), aminoAcid) != categorizer.bases.end()) ? true : false;
  int charge = isBase ? data.vectorOfChargesAndMultiplicities[fragmentIndex].first + 1
                      : data.vectorOfChargesAndMultiplicities[fragmentIndex].first - 1;

  Utils::FilesystemHelpers::createDirectories(
      Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(fragmentIndex), nonRefStateDir));

  std::string infoFilename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(fragmentIndex),
                                                                         nonRefStateDir, "info.dat");
  std::ofstream infoFile(infoFilename);
  infoFile << charge << "  " << data.vectorOfChargesAndMultiplicities[fragmentIndex].second << std::endl;
  infoFile.close();

  std::vector<int> superfluousHydrogens;
  std::list<int> neighbors;
  auto refStructure = *data.vectorOfStructures[fragmentIndex];
  if (data.vectorOfStructures.size() == 1)
    neighbors = data.listsOfNeighbors[criticalAtomIndex];
  else
    neighbors = data.listsOfNeighbors[fragmentIndex];

  for (auto& neighbor : neighbors) {
    if (data.fullStructure.at(neighbor).getElementType() == Utils::ElementType::H) {
      // get index of hydrogen atom in fragment
      int indexOfHInFragment = std::distance(data.atomIndexMapping.at(fragmentIndex).begin(),
                                             find(data.atomIndexMapping.at(fragmentIndex).begin(),
                                                  data.atomIndexMapping.at(fragmentIndex).end(), neighbor));
      superfluousHydrogens.push_back(indexOfHInFragment);
    }
  }
  TitrationHelper helper(settings);
  auto modifiedStructure =
      helper.changeProtonationState(refStructure, aminoAcid, isBase, criticalAtomIndex, superfluousHydrogens);
  Utils::ChemicalFileHandler::write(Utils::NativeFilenames::combinePathSegments(
                                        referenceDataDir, std::to_string(fragmentIndex), nonRefStateDir, "molecule.xyz"),
                                    modifiedStructure);

  std::vector<int> indicesOfAtomsInFragment;
  SwooseUtilities::FragmentationHelper::updateInformationForIndexMapping(refStructure, modifiedStructure,
                                                                         indicesOfAtomsInFragment);

  std::string constrFilename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(fragmentIndex),
                                                                           nonRefStateDir, "constr.dat");
  std::vector<int> constrainedAtoms;
  for (const auto& c : data.constrainedAtoms[fragmentIndex]) {
    constrainedAtoms.push_back(indicesOfAtomsInFragment.at(c));
  }
  writeConstrainedAtomsFile(constrainedAtoms, constrFilename);
}

void saveAdditionalStructuresForTitration(ParametrizationData& data, TitrationResults& results, int fragmentIndex,
                                          std::string referenceDataDir) {
  if (data.siteIspHSensitive.at(fragmentIndex)) {
    std::string filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(fragmentIndex),
                                                                       nonRefStateDir, "opt.xyz");
    try {
      auto optStructure = Utils::ChemicalFileHandler::read(filename).first;
      results.vectorOfOptimizedNonRefStructures.insert({fragmentIndex, std::make_unique<Utils::AtomCollection>(optStructure)});
    }
    catch (const std::exception& e) {
      results.vectorOfOptimizedNonRefStructures.insert({fragmentIndex, nullptr});
    }
  }
}

void parseElectronicEnergiesForTitration(ParametrizationData& data, TitrationResults& results, int fragmentIndex,
                                         std::string referenceDataDir, bool parseTurbomoleOutput,
                                         std::shared_ptr<Utils::Settings> settings) {
  // TODO: this is still very ugly
  if (data.siteIspHSensitive.at(fragmentIndex)) {
    // clean up
    results.electronicEnergies.erase(fragmentIndex);
    bool useThermochemistryData = settings->getBool(SwooseUtilities::SettingsNames::useThermoChemistryForTitration);
    try {
      double energyRef;
      double energyNonRef;
      if (parseTurbomoleOutput) {
        auto turbomoleOutputParser1 = getPreparedTurbomoleParser(referenceDataDir, fragmentIndex);
        Utils::ExternalQC::TurbomoleFiles outputFiles;
        Utils::ExternalQC::setCorrectTurbomoleFileNames(
            outputFiles,
            Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(fragmentIndex), nonRefStateDir));
        Utils::ExternalQC::TurbomoleMainOutputParser parser2(outputFiles);
        if (useThermochemistryData) {
          throw std::runtime_error("Turbomole calculator doesn't support thermochemistry. ");
        }
        energyRef = turbomoleOutputParser1.getEnergy();
        energyNonRef = parser2.getEnergy();
      }
      else {
        std::string filenameNonRef = Utils::NativeFilenames::combinePathSegments(
            referenceDataDir, std::to_string(fragmentIndex), nonRefStateDir, "orca.out");
        std::string filenameRef =
            Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(fragmentIndex), "orca.out");
        Utils::ExternalQC::OrcaMainOutputParser orcaMainOutputParser1(filenameRef);
        Utils::ExternalQC::OrcaMainOutputParser orcaMainOutputParser2(filenameNonRef);
        if (useThermochemistryData) {
          energyRef = orcaMainOutputParser1.getGibbsFreeEnergy();
          energyNonRef = orcaMainOutputParser2.getGibbsFreeEnergy();
        }
        else {
          energyRef = orcaMainOutputParser1.getEnergy();
          energyNonRef = orcaMainOutputParser2.getEnergy();
        }
      }
      if (results.electronicEnergies.count(fragmentIndex) == 1)
        throw std::runtime_error("Information for fragment " + std::to_string(fragmentIndex) + " already present.");

      results.electronicEnergies.insert({fragmentIndex, std::make_pair(energyRef, energyNonRef)});
    }
    catch (const std::exception& e) {
      throw std::runtime_error("No electronic energies provided for titration calculation of fragment " +
                               std::to_string(fragmentIndex) + " Error: " + e.what());
    }
  }
}

void readReferenceDataFromFiles(ParametrizationData& data, TitrationResults& titrationResults,
                                std::string referenceDataDir, std::shared_ptr<Utils::Settings> settings, Core::Log& log) {
  log.debug << "Reading optimized structures from disk..." << Core::Log::endl;
  data.vectorOfOptimizedStructures.resize(data.vectorOfStructures.size());
  bool readFromCsv = settings->getBool(SwooseUtilities::SettingsNames::useCsvInputFormat);
  auto referenceProgram = settings->getString(SwooseUtilities::SettingsNames::referenceProgram);
  bool parseTurbomoleOutput = false;

  if (!readFromCsv) {
    if (referenceProgram != SwooseUtilities::OptionNames::turbomoleOption &&
        referenceProgram != SwooseUtilities::OptionNames::orcaOption)
      throw std::runtime_error("Reading reference data from files is only available for Turbomole and Orca. For other "
                               "programs, please use the csv file format option.");
    if (referenceProgram == SwooseUtilities::OptionNames::turbomoleOption)
      parseTurbomoleOutput = true;
  }

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < int(data.vectorOfStructures.size()); ++i) {
    std::string filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "opt.xyz");
    try {
      auto optStructure = Utils::ChemicalFileHandler::read(filename).first;
      data.vectorOfOptimizedStructures[i] = std::make_unique<Utils::AtomCollection>(optStructure);
    }
    catch (const std::exception& e) {
      data.vectorOfOptimizedStructures[i] = nullptr;
    }
    if (settings->getBool(SwooseUtilities::SettingsNames::titrate)) {
      saveAdditionalStructuresForTitration(data, titrationResults, i, referenceDataDir);
      parseElectronicEnergiesForTitration(data, titrationResults, i, referenceDataDir, parseTurbomoleOutput, settings);
    }
  }

  log.debug << "Reading Hessians from disk..." << Core::Log::endl;
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < int(data.vectorOfStructures.size()); ++i) {
    // Check whether structure optimization succeeded
    if (data.vectorOfOptimizedStructures.at(i) == nullptr) {
      data.vectorOfHessians[i] = nullptr;
      continue;
    }

    std::string filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "hessian.csv");
    try {
      Utils::HessianMatrix hessianMatrix;
      if (readFromCsv)
        hessianMatrix = Utils::csvToMatrix(filename, ',');
      else if (parseTurbomoleOutput) {
        auto turbomoleOutputParser = getPreparedTurbomoleParser(referenceDataDir, i);
        hessianMatrix = turbomoleOutputParser.getHessian();
      }
      else {
        filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "orca.hess");
        hessianMatrix = Utils::ExternalQC::OrcaHessianOutputParser::getHessian(filename);
      }
      data.vectorOfHessians[i] = std::make_unique<Utils::HessianMatrix>(hessianMatrix);
    }
    catch (const std::exception& e) {
      data.vectorOfHessians[i] = nullptr;
    }
  }

  log.debug << "Reading atomic charges from disk..." << Core::Log::endl;
  data.atomicChargesForEachFragment.resize(data.vectorOfStructures.size());
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < int(data.vectorOfStructures.size()); ++i) {
    // Check whether structure optimization succeeded
    if (data.vectorOfOptimizedStructures.at(i) == nullptr) {
      data.atomicChargesForEachFragment[i].clear(); // Should already be empty, but still for consistency
      continue;
    }

    if (readFromCsv) {
      std::string filename =
          Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "atomic_charges.csv");
      try {
        Eigen::MatrixXd atomicChargesMatrix = Utils::csvToMatrix(filename, ',');
        std::vector<double> atomicCharges(
            atomicChargesMatrix.data(), atomicChargesMatrix.data() + atomicChargesMatrix.rows() * atomicChargesMatrix.cols());
        if (settings->getBool(SwooseUtilities::SettingsNames::convertChargesCm5))
          data.atomicChargesForEachFragment[i] =
              Utils::ChargeModel5::calculateCm5Charges(atomicCharges, *data.vectorOfOptimizedStructures.at(i));
        else
          data.atomicChargesForEachFragment[i] = atomicCharges;
      }
      catch (const std::exception& e) {
        data.atomicChargesForEachFragment[i].clear(); // Should already be empty, but still for consistency
      }
    }
    else if (settings->getBool(SwooseUtilities::SettingsNames::useGaussianOptionKey)) {
      std::string filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "cm5.out");
      try {
        Utils::ExternalQC::GaussianOutputParser gaussianParser(filename);
        data.atomicChargesForEachFragment[i] = gaussianParser.getCM5Charges();
      }
      catch (const std::exception& e) {
        data.atomicChargesForEachFragment[i].clear(); // Should already be empty, but still for consistency
      }
    }
    else {
      std::string filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "hirshfeld.out");
      try {
        // Note: atomic charges are Hirshfeld charges for Orca and Loewdin charges for Turbomole
        std::vector<double> charges;
        if (parseTurbomoleOutput) {
          auto turbomoleOutputParser = getPreparedTurbomoleParser(referenceDataDir, i);
          charges = turbomoleOutputParser.getLoewdinCharges();
        }
        else {
          Utils::ExternalQC::OrcaMainOutputParser orcaMainOutputParser(filename);
          charges = orcaMainOutputParser.getHirshfeldCharges();
        }
        data.atomicChargesForEachFragment[i] =
            Utils::ChargeModel5::calculateCm5Charges(charges, *data.vectorOfOptimizedStructures.at(i));
      }
      catch (const std::exception& e) {
        data.atomicChargesForEachFragment[i].clear(); // Should already be empty, but still for consistency
      }
    }
  }

  if (settings->getBool(SwooseUtilities::SettingsNames::refineConnectivity)) {
    log.debug << "Reading covalent bond orders from disk..." << Core::Log::endl;
    data.vectorOfBondOrderCollections.resize(data.vectorOfStructures.size()); // resize vector of bond orders
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < int(data.vectorOfStructures.size()); ++i) {
      std::string filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "bonds.out");
      if (readFromCsv)
        filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "connectivity.dat");
      try {
        if (readFromCsv) {
          auto listsOfNeighbors = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(filename);
          data.vectorOfBondOrderCollections[i] = std::make_unique<Utils::BondOrderCollection>(
              SwooseUtilities::TopologyUtils::generateBondOrderMatrixFromListsOfNeighbors(listsOfNeighbors));
        }
        else if (parseTurbomoleOutput) {
          auto turbomoleOutputParser = getPreparedTurbomoleParser(referenceDataDir, i);
          data.vectorOfBondOrderCollections[i] =
              std::make_unique<Utils::BondOrderCollection>(turbomoleOutputParser.getBondOrders());
        }
        else {
          Utils::ExternalQC::OrcaMainOutputParser orcaMainOutputParser(filename);
          data.vectorOfBondOrderCollections[i] =
              std::make_unique<Utils::BondOrderCollection>(orcaMainOutputParser.getBondOrders());
        }
      }
      catch (const std::exception& e) {
        data.vectorOfBondOrderCollections[i] = nullptr;
      }
    }
  }
}

Utils::ExternalQC::TurbomoleMainOutputParser getPreparedTurbomoleParser(const std::string& referenceDataDir, int fragmentIndex) {
  Utils::ExternalQC::TurbomoleFiles outputFiles;
  Utils::ExternalQC::setCorrectTurbomoleFileNames(
      outputFiles, Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(fragmentIndex)));
  Utils::ExternalQC::TurbomoleMainOutputParser parser(outputFiles);
  return parser;
}

void writeConstrainedAtomsFile(const std::vector<int>& constrainedAtoms, std::string& constrainedAtomsFile) {
  if (!constrainedAtoms.empty()) {
    std::string constraintsString;
    for (const auto& c : constrainedAtoms) {
      constraintsString += std::to_string(c);
      constraintsString += " ";
    }
    std::ofstream constrFile(constrainedAtomsFile);
    constrFile << constraintsString << std::endl;
    constrFile.close();
  }
}

} // namespace ReferenceCalculationsIO
} // namespace MMParametrization
} // namespace Scine
