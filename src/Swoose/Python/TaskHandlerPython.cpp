/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "../App/AppUtils/TaskManagement.h"
#include <Core/Log.h>
#include <Swoose/Utilities/ModuleLoader.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/Technical/UniqueIdentifier.h>
#include <pybind11/pybind11.h>
#include <boost/dll/runtime_symbol_info.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <vector>

using namespace Scine;

/*
 * Helper construct to make sure that the helper yaml file is always deleted even if exception is thrown.
 * Destructor takes care of deleting the temporary yaml file.
 */
struct YamlFileHandler {
  YamlFileHandler() {
    Utils::UniqueIdentifier filenameId;
    filenameIdString = filenameId.getStringRepresentation();
  }
  ~YamlFileHandler() {
    if (!yamlFilePath.empty())
      boost::filesystem::remove_all(yamlFilePath);
  }
  void dumpYamlFileAndUpdateFilename(const YAML::Node& yamlSettings) {
    yamlFilePath = "tmp_settings_" + filenameIdString + ".yaml";
    std::ofstream yamlFile(yamlFilePath);
    yamlFile << yamlSettings << std::endl;
    yamlFile.close();
  };
  std::string yamlFilePath;
  std::string filenameIdString;
};

void runTask(std::string mode, bool quantum, std::string structureFile, pybind11::kwargs kwargs) {
  // Get a module manager
  boost::filesystem::path thisFile = boost::dll::this_line_location();
  std::string packageDirectory = thisFile.parent_path().parent_path().string();
  setenv("SCINE_MODULE_PATH",
         (packageDirectory + "/scine_swoose:" + packageDirectory + "/scine_utilities:" + packageDirectory +
          "/scine_sparrow:" + packageDirectory + "/scine_xtb_wrapper")
             .c_str(),
         false);
  auto& manager = Core::ModuleManager::getInstance();
  SwooseUtilities::loadModules(manager, {"Swoose", "Orca", "Gaussian", "Turbomole"}, true); // essential modules
  SwooseUtilities::loadModules(manager, {"Sparrow", "Xtb"}, false);                         // optional modules

  bool verbose = false;
  bool hessianRequired = false;

  // Read settings from kwargs
  YAML::Node yamlSettings;
  for (auto item : kwargs) {
    std::string key = item.first.cast<std::string>();
    if (!key.compare("verbose")) {
      verbose = item.second.cast<bool>();
    }
    else if (!key.compare("hessian")) {
      hessianRequired = item.first.cast<bool>();
    }
    else if (pybind11::isinstance<pybind11::bool_>(item.second)) {
      bool value = item.second.cast<bool>();
      yamlSettings[key] = value;
    }
    else if (pybind11::isinstance<pybind11::str>(item.second)) {
      std::string value = item.second.cast<std::string>();
      yamlSettings[key] = value;
    }
    else if (pybind11::isinstance<pybind11::int_>(item.second)) {
      int value = item.second.cast<int>();
      yamlSettings[key] = value;
    }
    else if (pybind11::isinstance<pybind11::float_>(item.second)) {
      double value = item.second.cast<double>();
      yamlSettings[key] = value;
    }
    else if (pybind11::isinstance<pybind11::list>(item.second)) {
      std::vector<int> values;
      pybind11::list outputList = item.second.cast<pybind11::list>();
      for (auto item : outputList) {
        int value = item.cast<int>();
        values.push_back(value);
      }
      yamlSettings[key] = values;
    }
    else {
      throw std::runtime_error(key + " could not be converted from Python into C++ a variable, check its type!");
    }
  }

  // Create logger
  Core::Log log;
  if (verbose)
    log.debug.add("cout", Core::Log::coutSink());

  /*
   * Let the task manager run the tasks. For the direct mode in parametrization,
   * we need to dump the settings to a yaml file.
   */
  YamlFileHandler yamlFileHandler;
  if (mode == "parametrize" || mode == "select_qm") {
    yamlFileHandler.dumpYamlFileAndUpdateFilename(yamlSettings);
  }
  Swoose::TaskManagement::manageTasks(manager, mode, quantum, hessianRequired, structureFile, yamlSettings,
                                      yamlFileHandler.yamlFilePath, log);
}

void runParametrization(std::string structureFile, pybind11::kwargs kwargs) {
  runTask("parametrize", false, structureFile, kwargs);
}

void runMmCalculation(std::string structureFile, pybind11::kwargs kwargs) {
  runTask("calculate", false, structureFile, kwargs);
}

void runQmmmCalculation(std::string structureFile, pybind11::kwargs kwargs) {
  runTask("calculate", true, structureFile, kwargs);
}

void runMmOptimization(std::string structureFile, pybind11::kwargs kwargs) {
  runTask("optimize", false, structureFile, kwargs);
}

void runQmmmOptimization(std::string structureFile, pybind11::kwargs kwargs) {
  runTask("optimize", true, structureFile, kwargs);
}

void runMdSimulation(std::string structureFile, pybind11::kwargs kwargs) {
  runTask("md", false, structureFile, kwargs);
}

void runMdSimulationWithQmmm(std::string structureFile, pybind11::kwargs kwargs) {
  runTask("md", true, structureFile, kwargs);
}

void runQmRegionSelection(std::string structureFile, pybind11::kwargs kwargs) {
  runTask("select_qm", false, structureFile, kwargs);
}

void init_tasks(pybind11::module& m) {
  m.def("parametrize", &runParametrization, pybind11::arg("structure file"),
        "Parametrizes the SFAM model. Settings can be set as keyword arguments.");
  m.def("calculate_mm", &runMmCalculation, pybind11::arg("structure file"),
        "Calculation with a molecular mechanics model. Settings can be set as keyword arguments.");
  m.def("calculate_qmmm", &runQmmmCalculation, pybind11::arg("structure file"),
        "Calculation with the QM/MM hybrid model. Settings can be set as keyword arguments.");
  m.def("optimize_mm", &runMmOptimization, pybind11::arg("structure file"),
        "Structure optimization with a molecular mechanics model. Settings can be set as keyword arguments.");
  m.def("optimize_qmmm", &runQmmmOptimization, pybind11::arg("structure file"),
        "Structure optimization with the QM/MM hybrid model. Settings can be set as keyword arguments.");
  m.def("md_simulate_mm", &runMdSimulation, pybind11::arg("structure file"),
        "Molecular dynamics simulation with a molecular mechanics model. Settings can be set as keyword "
        "arguments.");
  m.def("md_simulate_qmmm", &runMdSimulationWithQmmm, pybind11::arg("structure file"),
        "Molecular dynamics simulation with a QM/MM model. Settings can be set as keyword "
        "arguments.");
  m.def("select_qm_region", &runQmRegionSelection, pybind11::arg("structure file"),
        "Automated QM region selection. Settings can be set as keyword arguments.");
}
