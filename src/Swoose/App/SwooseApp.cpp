/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AppUtils/CommandLineOptions.h"
#include "AppUtils/TaskManagement.h"
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Swoose/Utilities/ModuleLoader.h>
#include <Utils/GeometryOptimization/QmmmGeometryOptimizer.h>
#include <Utils/IO/Yaml.h>
#include <Utils/MolecularDynamics/MolecularDynamics.h>
#include <yaml-cpp/yaml.h>
#include <chrono>

using namespace Scine;

// Declaration of two printing functions
void printHeader(Core::Log& log);
void printTiming(Core::Log& log, int start, int end);

int main(int argc, char* argv[]) {
  // Create logger
  Core::Log log;
  // Header
  printHeader(log);

  Swoose::CommandLineOptions commandLineParser(argc, argv);

  if (commandLineParser.helpRequired() || argc == 1) {
    commandLineParser.printHelp(log);
    return 0;
  }

  if (commandLineParser.debugLoggingRequired()) {
    log.debug.add("cout", Core::Log::coutSink());
  }

  // Get settings in YAML format and molecular structure file
  auto settingsFilePath = commandLineParser.getSettingsFile();
  YAML::Node yamlSettings;

  if (!settingsFilePath.empty())
    yamlSettings = YAML::LoadFile(settingsFilePath);
  else
    log.output << "No settings file was specified. Working with the defaults...\n" << Core::Log::endl;

  auto structureFile = commandLineParser.getStructureFile();

  // Get the information whether QM/MM is the method of choice
  bool quantum = commandLineParser.quantumCalculationRequired();

  // Get a module manager
  auto& manager = Core::ModuleManager::getInstance();
  SwooseUtilities::loadModules(manager, {"Swoose", "Orca", "Gaussian"}, true); // essential modules
  SwooseUtilities::loadModules(manager, {"Sparrow", "Xtb"}, false);            // optional modules

  // Handle the different possible modes of the app
  auto mode = commandLineParser.getMode();

  // Hessian calculation required?
  auto hessianRequired = commandLineParser.hessianRequired();

  auto start =
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  // Let the task manager run the tasks.
  Swoose::TaskManagement::manageTasks(manager, mode, quantum, hessianRequired, structureFile, yamlSettings,
                                      settingsFilePath, log);

  auto end =
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  // Timing
  printTiming(log, start, end);
  return 0;
}

void printHeader(Core::Log& log) {
  log.output << R"(#=========================================#)" << Core::Log::nl;
  log.output << R"(|    ____                                 |)" << Core::Log::nl;
  log.output << R"(|   / ___|_      _____   ___  ___  ___    |)" << Core::Log::nl;
  log.output << R"(|   \___ \ \ /\ / / _ \ / _ \/ __|/ _ \   |)" << Core::Log::nl;
  log.output << R"(|    ___) \ V  V / (_) | (_) \__ \  __/   |)" << Core::Log::nl;
  log.output << R"(|   |____/ \_/\_/ \___/ \___/|___/\___|   |)" << Core::Log::nl;
  log.output << R"(|                                         |)" << Core::Log::nl;
  log.output << R"(#=========================================#)" << Core::Log::nl;
  log.output << Core::Log::endl;
}

void printTiming(Core::Log& log, int start, int end) {
  int duration = end - start;
  int hour = 0;
  int min = 0;
  int sec = 0;
  int millisec = 0;

  hour = duration / 3600000;
  duration = duration % 3600000;
  min = duration / 60000;
  duration = duration % 60000;
  sec = duration / 1000;
  duration = duration % 1000;
  millisec = duration;

  log.output << "Total job duration:" << Core::Log::nl << hour << " hours, " << min << " minutes, " << sec
             << " seconds, " << millisec << " milliseconds." << Core::Log::endl;
}
