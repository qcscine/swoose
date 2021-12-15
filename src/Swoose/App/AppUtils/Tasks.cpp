/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Tasks.h"
#include <Core/Log.h>
#include <Swoose/QMMM/QmRegionSelection/QmRegionSelector.h>
#include <Swoose/Utilities/SettingsNames.h>
#include <Utils/GeometryOptimization/QmmmGeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/MolecularDynamics/MolecularDynamics.h>
#include <Utils/Optimizer/GradientBased/Bfgs.h>
#include <Utils/Optimizer/GradientBased/GradientBasedCheck.h>
#include <boost/filesystem.hpp>

namespace Scine {
namespace Swoose {
namespace Tasks {

namespace {
constexpr const char* defaultOptimizationResultsDirectory = "opt_results";
constexpr const char* optimizationDirectoryKey = "results_directory";
constexpr const char* optimizationStructureFile = "opt_structure.xyz";
constexpr const char* optimizationTrajectoryFile = "opt_trajectory.xyz";
constexpr const char* hessianCsvFilename = "hessian.csv";
} // namespace

void runMMCalculationTask(Core::Calculator& calculator, const std::string& structureFile,
                          Utils::PropertyList properties, Core::Log& log) {
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(structureFile).first;

  calculator.setStructure(structure);
  calculator.setRequiredProperties(properties);

  log.output << "Starting an MM calculation for the molecule: " << structureFile << Core::Log::nl;
  log.output << "MM model: " << calculator.name() << Core::Log::endl;
  const Utils::Results& results = calculator.calculate("MM calculation");
  double energy = results.get<Utils::Property::Energy>();
  Utils::GradientCollection gradients = results.get<Utils::Property::Gradients>() * Utils::Constants::kCalPerMol_per_hartree;

  log.output << "Calculation done. Results: " << Core::Log::endl;
  log.output << Core::Log::nl << "# Energy (kcal/mol): " << energy * Utils::Constants::kCalPerMol_per_hartree
             << Core::Log::endl;
  log.output << Core::Log::nl << "# Gradients (kcal/(mol*bohr)):" << Core::Log::endl;
  log.output << [&gradients](std::ostream& os) { Utils::matrixPrettyPrint(os, gradients); };
  log.output << Core::Log::endl;

  if (properties.containsSubSet(Utils::Property::Hessian)) {
    Utils::HessianMatrix hessian = results.get<Utils::Property::Hessian>() * Utils::Constants::kCalPerMol_per_hartree;
    std::ofstream hessianFile(hessianCsvFilename);
    Utils::matrixToCsv(hessianFile, hessian, ',');
    hessianFile.close();
    log.output << Core::Log::nl << "# Hessian (kcal/(mol*bohr^2)) has been written to file: " << hessianCsvFilename
               << Core::Log::nl << Core::Log::endl;
  }
}

void runQmmmCalculationTask(Core::Calculator& calculator, const std::string& structureFile,
                            Utils::PropertyList properties, Core::Log& log, const YAML::Node& yamlNode) {
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(structureFile).first;

  // Set settings
  Utils::nodeToSettings(calculator.settings(), yamlNode, false);
  // Set structure
  calculator.setStructure(structure);
  // Set required properties
  calculator.setRequiredProperties(properties);

  log.output << "Starting a QM/MM calculation for the molecule: " << structureFile << Core::Log::endl;
  const Utils::Results& results = calculator.calculate("QM/MM calculation");
  double energy = results.get<Utils::Property::Energy>();
  Utils::GradientCollection gradients = results.get<Utils::Property::Gradients>() * Utils::Constants::kCalPerMol_per_hartree;

  log.output << "Calculation done. Results: " << Core::Log::endl;
  log.output << Core::Log::nl << "# Energy (Hartree): " << std::setprecision(10) << energy << Core::Log::endl;
  log.output << Core::Log::nl << "# Gradients (kcal/(mol*bohr)):" << std::setprecision(5) << Core::Log::endl;
  log.output << [&gradients](std::ostream& os) { Utils::matrixPrettyPrint(os, gradients); };
  log.output << Core::Log::endl;
}

void runSFAMParametrizationTask(Core::MMParametrizer& parametrizer, const std::string& structureFile, Core::Log& log) {
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(structureFile).first;
  parametrizer.parametrize(structure);
  log.output << Core::Log::nl << "Parametrization done." << Core::Log::endl;
}

void runMDSimulationTask(Utils::MolecularDynamics& md, const std::string& structureFile, Core::Log& log) {
  Utils::AtomCollection startStructure = Utils::ChemicalFileHandler::read(structureFile).first;
  log.output << "Starting MD simulation..." << Core::Log::endl;
  try {
    md.performMDSimulation(startStructure, log);
  }
  catch (const std::exception& e) {
    const std::string failedTrajectoryFilename = "MD_trajectory_failed.xyz";
    auto failedTrajectory = md.getMolecularTrajectory();
    Utils::MolecularTrajectoryIO::write(Utils::MolecularTrajectoryIO::format::xyz, failedTrajectoryFilename, failedTrajectory);
    throw std::runtime_error(std::string(e.what()) + "\n" +
                             "Trajectory of failed MD has been written to: " + std::string(failedTrajectoryFilename));
  }

  log.output << "All MD steps have been successfully completed." << Core::Log::endl;

  auto trajectory = md.getMolecularTrajectory();
  const std::string trajectoryFilename = "MD_trajectory.xyz";
  Utils::MolecularTrajectoryIO::write(Utils::MolecularTrajectoryIO::format::xyz, trajectoryFilename, trajectory);
  log.output << "Trajectory has been written to the following file: " << trajectoryFilename << Core::Log::endl;

  auto energies = trajectory.getEnergies();
  const std::string energiesFilename = "MD_energies.dat";
  std::ofstream energiesFile(energiesFilename);
  energiesFile << "Energies for the MD snapshots (kcal/mol):" << Core::Log::endl;
  for (const auto& energy : energies) {
    energiesFile << energy * Utils::Constants::kCalPerMol_per_hartree << Core::Log::endl;
  }
  log.output << "Energies have been written to the following file: " << energiesFilename << Core::Log::endl;
}

void runMMOptimizationTask(Core::Calculator& calculator, Utils::GeometryOptimizerBase& optimizer,
                           const std::string& structureFile, Core::Log& log, const YAML::Node& yamlNode) {
  std::string resultsDir = defaultOptimizationResultsDirectory;
  if (yamlNode[optimizationDirectoryKey])
    resultsDir = yamlNode[optimizationDirectoryKey].as<std::string>();

  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(structureFile).first;

  Utils::FilesystemHelpers::createDirectories(resultsDir);

  // Trajectory stream
  boost::filesystem::remove(Utils::NativeFilenames::combinePathSegments(resultsDir, optimizationTrajectoryFile));
  std::ofstream trajectory(Utils::NativeFilenames::combinePathSegments(resultsDir, optimizationTrajectoryFile),
                           std::ofstream::out);
  Utils::XyzStreamHandler writer;
  // Variable for the observer
  double oldEnergy = 0.0;

  // Observer function
  auto obsFunc = [&](const int& cycle, const double& energy, const Eigen::VectorXd& /* params */) {
    if (cycle == 1) {
      log.output.printf("%s\n", "Starting optimization cycles...");
      log.output.printf("%7s %16s %16s\n", "Cycle", "Energy", "Energy Diff.");
    }
    log.output.printf("%7d %+16.9f %+16.9f\n", cycle, energy, energy - oldEnergy);
    oldEnergy = energy;
    auto currStructure = calculator.getStructure();
    writer.write(trajectory, *currStructure);
  };

  // Add observer
  optimizer.addObserver(obsFunc);

  auto cycles = optimizer.optimize(structure, log);

  trajectory.close();
  Utils::ChemicalFileHandler::write(Utils::NativeFilenames::combinePathSegments(resultsDir, optimizationStructureFile),
                                    structure);

  auto maxAllowedCycles = optimizer.getSettings().getInt(Utils::GradientBasedCheck::gconvMaxIterKey);
  if (cycles == maxAllowedCycles)
    throw std::runtime_error("The structure optimization did not converge within the maximum number of cycles. This "
                             "number can be increased with the setting: " +
                             std::string(Utils::GradientBasedCheck::gconvMaxIterKey));

  log.output << "MM structure optimization done. Results are written to directory: " << resultsDir << Core::Log::endl;
}

template<class OptimizerType>
void runQmmmOptimizationTask(Core::Calculator& calculator, Utils::QmmmGeometryOptimizer<OptimizerType>& optimizer,
                             const std::string& structureFile, Core::Log& log, const YAML::Node& yamlNode) {
  std::string resultsDir = defaultOptimizationResultsDirectory;
  if (yamlNode[optimizationDirectoryKey])
    resultsDir = yamlNode[optimizationDirectoryKey].as<std::string>();

  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(structureFile).first;

  // Set settings
  Utils::nodeToSettings(calculator.settings(), yamlNode, true);
  // Set structure
  calculator.setStructure(structure);

  Utils::FilesystemHelpers::createDirectories(resultsDir);

  // Trajectory stream
  boost::filesystem::remove(Utils::NativeFilenames::combinePathSegments(resultsDir, optimizationTrajectoryFile));
  std::ofstream trajectory(Utils::NativeFilenames::combinePathSegments(resultsDir, optimizationTrajectoryFile),
                           std::ofstream::out);
  Utils::XyzStreamHandler writer;
  // Variable for the observer
  double oldEnergy = 0.0;

  auto obsFunc = [&](const int& cycle, const double& energy, const Eigen::VectorXd& /* params */) {
    if (cycle == 1) {
      log.output.printf("%s\n", "Starting full system optimization cycles...");
      log.output.printf("%7s %16s %16s\n", "Cycle", "Energy", "Energy Diff.");
    }
    log.output.printf("%7d %+16.9f %+16.9f\n", cycle, energy, energy - oldEnergy);
    oldEnergy = energy;
    auto currStructure = calculator.getStructure();
    writer.write(trajectory, *currStructure);
  };

  // Add observer
  optimizer.addObserver(obsFunc);

  auto cycles = optimizer.optimize(structure, log);

  trajectory.close();
  Utils::ChemicalFileHandler::write(Utils::NativeFilenames::combinePathSegments(resultsDir, optimizationStructureFile),
                                    structure);

  auto maxAllowedCycles =
      optimizer.getSettings().getInt(Utils::QmmmGeometryOptimizer<Utils::Bfgs>::qmmmOptMaxMacroiterationsKey) *
      optimizer.getSettings().getInt(Utils::QmmmGeometryOptimizer<Utils::Bfgs>::qmmmOptMaxFullMicroiterationsKey);
  if (cycles == maxAllowedCycles)
    throw std::runtime_error("The structure optimization did not converge within the maximum number of cycles.");

  log.output << "QM/MM structure optimization done. Results are written to directory: " << resultsDir << Core::Log::endl;
}

void runQmRegionSelectionTask(const std::string& structureFile, Core::Log& log, YAML::Node& yamlNode,
                              std::string yamlSettingsPath) {
  log.output << "Starting automated selection of QM region..." << Core::Log::nl << Core::Log::endl;
  Qmmm::QmRegionSelector qmRegionSelector;
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(structureFile).first;
  qmRegionSelector.setLog(log);
  Utils::nodeToSettings(qmRegionSelector.settings(), yamlNode, true);
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::yamlSettingsFilePath, yamlSettingsPath);

  // Generate QM region
  qmRegionSelector.generateQmRegion(structure);
  auto result = qmRegionSelector.getQmRegionIndices();
  log.output << "Indices of QM atoms: " << Core::Log::endl;
  log.output << "[";
  for (int i = 0; i < result.size(); ++i) {
    if (result[i] >= 0) {
      log.output << result[i];
      if (i != result.size() - 1)
        log.output << ", ";
    }
  }
  log.output << "]" << Core::Log::nl << Core::Log::endl;
  auto optimalQmRegion = qmRegionSelector.getQmRegionStructure();
  auto chargeAndMultiplicity = qmRegionSelector.getQmRegionChargeAndMultiplicity();
  log.output << "Molecular charge of QM region: " << chargeAndMultiplicity.first << Core::Log::nl;
  log.output << "Spin multiplicity of QM region: " << chargeAndMultiplicity.second << Core::Log::endl;
  log.output << "Optimal QM region is written to: optimal_qm_region.xyz" << Core::Log::endl;
  Utils::ChemicalFileHandler::write("optimal_qm_region.xyz", optimalQmRegion);
}

// Explicit instantiation such that the function implementation can stay in the .cpp file.
template void
runQmmmOptimizationTask<Utils::Bfgs>(Core::Calculator& calculator, Utils::QmmmGeometryOptimizer<Utils::Bfgs>& optimizer,
                                     const std::string& structureFile, Core::Log& log, const YAML::Node& yamlNode);

} // namespace Tasks
} // namespace Swoose
} // namespace Scine