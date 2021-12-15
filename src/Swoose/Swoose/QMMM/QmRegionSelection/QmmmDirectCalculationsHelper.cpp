/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmmmDirectCalculationsHelper.h"
#include "QmRegionSelector.h"
#include "QmRegionSelectorSettings.h"
#include "QmmmReferenceDataManager.h"
#include "SymmetryScores.h"
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Swoose/QMMM/QmmmCalculator.h>
#include <Swoose/Utilities/CalculatorOptions.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/IO/Yaml.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <yaml-cpp/yaml.h>
#include <thread>

namespace Scine {
namespace Qmmm {

QmmmDirectCalculationsHelper::QmmmDirectCalculationsHelper(const Utils::Settings& settings, Core::Log& log,
                                                           const Utils::AtomCollection& structure,
                                                           const std::vector<QmmmModel>& qmmmModelCandidates,
                                                           const std::vector<QmmmModel>& qmmmReferenceModels,
                                                           const QmmmData& qmmmData)
  : settings_(settings),
    log_(log),
    structure_(structure),
    qmmmModelCandidates_(qmmmModelCandidates),
    qmmmReferenceModels_(qmmmReferenceModels),
    qmmmData_(qmmmData) {
}

std::vector<ForcesCollection> QmmmDirectCalculationsHelper::calculateForces() {
  // Initialize some variables
  auto numCandidateModels = qmmmModelCandidates_.size();
  auto numModels = numCandidateModels + qmmmReferenceModels_.size();
  std::vector<ForcesCollection> forces(numModels);

  // Get yaml settings for QM/MM
  auto yamlSettingsPath = settings_.getString(SwooseUtilities::SettingsNames::yamlSettingsFilePath);
  YAML::Node yamlNodeQmmmSettings;
  if (!yamlSettingsPath.empty())
    yamlNodeQmmmSettings = YAML::LoadFile(yamlSettingsPath);

  // Get the models
  const auto qmCalculatorModelAndModule = SwooseUtilities::getChosenQmCalculatorOption(yamlNodeQmmmSettings);
  const auto mmModel = SwooseUtilities::getChosenMMCalculatorOption(yamlNodeQmmmSettings);

  // Get a module manager
  auto& manager = Core::ModuleManager::getInstance();
  // See whether Swoose has to be loaded as a module and whether the MM calculator is already available.
  try {
    auto mmCalculatorFirstTest = manager.get<Core::Calculator>(mmModel);
  }
  catch (const std::exception& e) {
    if (!manager.moduleLoaded("Swoose"))
      manager.load("swoose");
    try {
      auto mmCalculatorSecondTest = manager.get<Core::Calculator>(mmModel);
    }
    catch (const std::exception& e) {
      throw std::runtime_error("MM calculator could not be loaded via the module system.");
    }
  }

  // Get max allowed symmetry score
  const double maxAllowedSymmetryScore = calculateMaxAllowedSymmetryScore(qmmmData_, settings_);

#pragma omp parallel for
  for (int i = 0; i < numModels; ++i) {
    int molecularCharge, spinMultiplicity;
    if (i < numCandidateModels) {
      molecularCharge = qmmmModelCandidates_.at(i).molecularCharge;
      spinMultiplicity = qmmmModelCandidates_.at(i).spinMultiplicity;
    }
    else {
      molecularCharge = qmmmReferenceModels_.at(i - numCandidateModels).molecularCharge;
      spinMultiplicity = qmmmReferenceModels_.at(i - numCandidateModels).spinMultiplicity;
    }

    if (qmmmData_.symmetryScores.at(i) > maxAllowedSymmetryScore) {
#pragma omp critical
      log_.debug << "Reference calculation for model " + std::to_string(i) + " skipped due to symmetry score."
                 << Core::Log::endl;
      continue;
    }

    try {
      // Get MM and QM calculators
      std::shared_ptr<Core::Calculator> qmCalculator, mmCalculator;
      try {
        qmCalculator = manager.get<Core::Calculator>(qmCalculatorModelAndModule.first, qmCalculatorModelAndModule.second);
        mmCalculator = manager.get<Core::Calculator>(mmModel);
      }
      catch (const std::runtime_error& e) {
        throw std::runtime_error(
            "The QM or MM calculator could not be loaded via the module system.\nCheck: (i) whether the requested "
            "calculators are available (see manual), (ii) that you have installed all "
            "relevant modules, and (iii) that all necessary environment variables are set (e.g., ORCA_BINARY_PATH or "
            "TURBODIR).");
      }

      // Create and configure QM/MM calculator
      Qmmm::QmmmCalculator calculator;

      // Set logger
      Core::Log warningLog = Core::Log::silent();
      warningLog.warning.add("cerr", Core::Log::cerrSink());
      warningLog.error.add("cerr", Core::Log::cerrSink());
      calculator.setLog(warningLog);

      calculator.setUnderlyingCalculators(qmCalculator, mmCalculator);
      Utils::PropertyList properties = Utils::Property::Energy | Utils::Property::Gradients;

      // Forward QM/MM settings to calculator
      Utils::nodeToSettings(calculator.settings(), yamlNodeQmmmSettings, true);
      calculator.settings().modifyInt(Utils::SettingsNames::spinMultiplicity, spinMultiplicity);
      calculator.settings().modifyInt(Utils::SettingsNames::molecularCharge, molecularCharge);
      calculator.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath,
                                         settings_.getString(SwooseUtilities::SettingsNames::connectivityFilePath));
      calculator.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath,
                                         settings_.getString(SwooseUtilities::SettingsNames::parameterFilePath));

      // Add qm_atoms setting
      std::vector<int> qmAtomIndices;
      auto& qmAtomIndicesRef = i < numCandidateModels ? qmmmModelCandidates_.at(i).qmAtomIndices
                                                      : qmmmReferenceModels_.at(i - numCandidateModels).qmAtomIndices;
      for (int idx : qmAtomIndicesRef) {
        if (idx >= 0) // Indices of -1 indicate a hydrogen link atom internally and are therefore ignored.
          qmAtomIndices.push_back(idx);
      }
      calculator.settings().modifyIntList(SwooseUtilities::SettingsNames::qmAtomsList, qmAtomIndices);

      calculator.setStructure(structure_);
      calculator.setRequiredProperties(properties);
      const Utils::Results& results = calculator.calculate("QM/MM calculation");
      Utils::GradientCollection gradients = results.get<Utils::Property::Gradients>();
      forces.at(i) = -gradients;
    }
    catch (const std::exception& e) {
#pragma omp critical
      {
        log_.output << "Reference calculation for model " + std::to_string(i) + " failed.\nError Message:\n"
                    << e.what() << Core::Log::nl << Core::Log::endl;
      }
    }
  }
  return forces;
}

} // namespace Qmmm
} // namespace Scine
