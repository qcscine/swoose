/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ReferenceCalculationsIO.h"
#include "../ParametrizationData.h"
#include <Core/Log.h>
#include <Swoose/MMParametrization/MMParametrizationSettings.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/ExternalQC/Gaussian/GaussianOutputParser.h>
#include <Utils/ExternalQC/Orca/OrcaHessianOutputParser.h>
#include <Utils/ExternalQC/Orca/OrcaMainOutputParser.h>
#include <Utils/ExternalQC/Turbomole/TurbomoleMainOutputParser.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Properties/AtomicCharges/ChargeModel5.h>
#include <fstream>

namespace Scine {
namespace MMParametrization {
namespace ReferenceCalculationsIO {

void writeXyzFiles(ParametrizationData& data, std::string referenceDataDir) {
  // Write the structures to disk as xyz files along with information about charge, multiplicity and constrained atoms.
  for (int i = 0; i < data.vectorOfStructures.size(); ++i) {
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
    if (!constrainedAtoms.empty()) {
      std::string constraintsString;
      for (const auto& c : constrainedAtoms) {
        constraintsString += std::to_string(c);
        constraintsString += " ";
      }
      std::string constrFilename =
          Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "constr.dat");
      std::ofstream constrFile(constrFilename);
      constrFile << constraintsString << std::endl;
      constrFile.close();
    }
  }
}

void readReferenceDataFromFiles(ParametrizationData& data, std::string referenceDataDir,
                                std::shared_ptr<Utils::Settings> settings, Core::Log& log) {
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
  for (int i = 0; i < data.vectorOfStructures.size(); ++i) {
    std::string filename = Utils::NativeFilenames::combinePathSegments(referenceDataDir, std::to_string(i), "opt.xyz");
    try {
      auto optStructure = Utils::ChemicalFileHandler::read(filename).first;
      data.vectorOfOptimizedStructures[i] = std::make_unique<Utils::AtomCollection>(optStructure);
    }
    catch (const std::exception& e) {
      data.vectorOfOptimizedStructures[i] = nullptr;
    }
  }

  log.debug << "Reading Hessians from disk..." << Core::Log::endl;
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < data.vectorOfStructures.size(); ++i) {
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
  for (int i = 0; i < data.vectorOfStructures.size(); ++i) {
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
    for (int i = 0; i < data.vectorOfStructures.size(); ++i) {
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

} // namespace ReferenceCalculationsIO
} // namespace MMParametrization
} // namespace Scine