/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
#include <utility>

namespace Scine {
namespace Qmmm {

QmmmDirectCalculationsHelper::QmmmDirectCalculationsHelper(std::shared_ptr<QmmmCalculator> qmmmCalculator,
                                                           const Utils::Settings& settings, Core::Log& log,
                                                           const Utils::AtomCollection& structure,
                                                           const std::vector<QmmmModel>& qmmmModelCandidates,
                                                           const std::vector<QmmmModel>& qmmmReferenceModels,
                                                           const QmmmData& qmmmData)
  : settings_(settings),
    log_(log),
    structure_(structure),
    qmmmModelCandidates_(qmmmModelCandidates),
    qmmmReferenceModels_(qmmmReferenceModels),
    qmmmData_(qmmmData),
    qmmmCalculator_(std::move(qmmmCalculator)) {
}

std::vector<ForcesCollection> QmmmDirectCalculationsHelper::calculateForces() {
  // Initialize some variables
  int numCandidateModels = qmmmModelCandidates_.size();
  int numModels = numCandidateModels + qmmmReferenceModels_.size();
  std::vector<ForcesCollection> forces(numModels);

  // Get max allowed symmetry score
  const double maxAllowedSymmetryScore = calculateMaxAllowedSymmetryScore(qmmmData_, settings_);

  // setup cloned calculators outside of OMP region, because they might not be thread safe
  std::vector<std::shared_ptr<QmmmCalculator>> calculators;
  auto nThreads = omp_get_max_threads();
  calculators.reserve(nThreads);
  for (int i = 0; i < nThreads; ++i) {
    calculators.push_back(qmmmCalculator_->clone());
  }

#pragma omp parallel for
  for (long unsigned int i = 0; i < numModels; ++i) {
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

    auto& calculator = calculators.at(omp_get_thread_num());

    try {
      // Set logger
      Utils::CalculationRoutines::setLog(*calculator, true, true, false);

      const Utils::PropertyList properties = Utils::Property::Energy | Utils::Property::Gradients;

      // Forward QM/MM settings to calculator
      calculator->settings().modifyInt(Utils::SettingsNames::spinMultiplicity, spinMultiplicity);
      calculator->settings().modifyInt(Utils::SettingsNames::molecularCharge, molecularCharge);
      calculator->settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath,
                                          settings_.getString(SwooseUtilities::SettingsNames::connectivityFilePath));
      calculator->settings().modifyString(Utils::SettingsNames::parameterFilePath,
                                          settings_.getString(Utils::SettingsNames::parameterFilePath));

      // Add qm_atoms setting
      std::vector<int> qmAtomIndices;
      const auto& qmAtomIndicesRef = i < numCandidateModels ? qmmmModelCandidates_.at(i).qmAtomIndices
                                                            : qmmmReferenceModels_.at(i - numCandidateModels).qmAtomIndices;
      for (const int idx : qmAtomIndicesRef) {
        if (idx >= 0) // Indices of -1 indicate a hydrogen link atom internally and are therefore ignored.
          qmAtomIndices.push_back(idx);
      }
      calculator->settings().modifyIntList(Utils::SettingsNames::qmAtomsList, qmAtomIndices);
      calculator->setStructure(structure_);
      calculator->setRequiredProperties(properties);
      const Utils::Results& results = calculator->calculate("QM/MM calculation");
      const Utils::GradientCollection gradients = results.get<Utils::Property::Gradients>();
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
