/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmRegionCandidateGenerator.h"
#include "QmRegionSelector.h"
#include "QmRegionSelectorSettings.h"
#include <Swoose/Utilities/AtomicInformationReader.h>
#include <Swoose/Utilities/FragmentAnalyzer.h>
#include <Swoose/Utilities/SubsystemGenerator.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <limits>

namespace Scine {
namespace Qmmm {

void QmRegionCandidateGenerator::generateQmRegionCandidates(std::vector<QmmmModel>& qmmmModelCandidates,
                                                            std::vector<QmmmModel>& qmmmReferenceModels,
                                                            const Utils::AtomCollection& fullStructure,
                                                            const Utils::BondOrderCollection& bondOrders,
                                                            const Utils::Settings& settings, Core::Log& log) {
  // Get from settings
  using namespace SwooseUtilities::SettingsNames;
  std::string atomicInfoFilePath = settings.getString(atomicInformationFile);
  int numStructures = settings.getInt(numAttemptsPerRadius);
  double probability = settings.getDouble(cuttingProbability);
  double initialRadius = settings.getDouble(initialRadiusForQmRegionSelection);
  int centerAtom = settings.getInt(qmRegionCenterAtom);
  int minSize = settings.getInt(qmRegionCandidateMinSize);
  int maxSize = settings.getInt(qmRegionCandidateMaxSize);
  int maxSizeRef = settings.getInt(qmRegionRefMaxSize);
  int minSizeRef = static_cast<int>(0.95 * maxSizeRef); // TODO: Also setting?
  int randomSeed = settings.getInt(qmRegionSelectionRandomSeed);

  if (maxSize < minSize)
    throw std::runtime_error("The maximum QM region size is smaller than the minimum one.");
  if (maxSizeRef < maxSize)
    throw std::runtime_error("The maximum QM region size is larger than the maximum size of the reference QM regions.");

  // If the full structure can be the reference, it should be the only reference.
  if (maxSizeRef >= fullStructure.size()) {
    maxSizeRef = fullStructure.size();
    minSizeRef = fullStructure.size();
  }

  std::map<int, int> formalCharges;
  std::map<int, int> unpairedElectrons;

  if (!atomicInfoFilePath.empty()) {
    SwooseUtilities::AtomicInformationReader reader(log);
    reader.read(atomicInfoFilePath, formalCharges, unpairedElectrons, fullStructure.size());
  }

  SwooseUtilities::FragmentAnalyzer fragmentAnalyzer(formalCharges, unpairedElectrons);
  SwooseUtilities::SubsystemGenerator generator(fullStructure, bondOrders, fragmentAnalyzer, 0.5,
                                                std::numeric_limits<int>::max(), log, randomSeed, probability);

  if (probability == 1.0) {
    log.output << "Because a cutting probability of 1.0 was selected, only one QM region will be generated based on "
                  "the provided initial radius setting.\nIf the resulting QM region is too large for your purposes, "
                  "re-try with a smaller value for the initial radius."
               << Core::Log::nl << Core::Log::endl;
    generateSingleQmRegionNonRandomly(qmmmModelCandidates, generator, fragmentAnalyzer, initialRadius, centerAtom);
    return;
  }

  std::vector<std::vector<int>> allQmAtomIndices;
  double r = initialRadius;
  int numExceededRefSizeLimit = 0;
  while (numExceededRefSizeLimit != numStructures) { // Stop if every qm region for given r is too large
    numExceededRefSizeLimit = 0;
    for (int k = 0; k < numStructures; ++k) {
      std::vector<int> indices;
      Utils::AtomCollection qmRegion =
          generator.generateSubsystem(centerAtom, indices, r * Utils::Constants::bohr_per_angstrom);
      std::sort(indices.begin(), indices.end());

      if (qmRegion.size() == fullStructure.size() || qmRegion.size() > maxSizeRef) // TODO: Add buffer, e.g., 20%?
        numExceededRefSizeLimit++;

      bool hasAcceptableSize = qmRegion.size() >= minSize && qmRegion.size() <= maxSize;
      bool isReferenceModel = qmRegion.size() >= minSizeRef && qmRegion.size() <= maxSizeRef;
      if (!hasAcceptableSize && !isReferenceModel)
        continue;
      bool unique = std::find(allQmAtomIndices.begin(), allQmAtomIndices.end(), indices) == allQmAtomIndices.end();
      if (!unique)
        continue;
      auto molecularCharge = fragmentAnalyzer.getMolecularCharge();
      auto spinMultiplicity = fragmentAnalyzer.getSpinMultiplicity();
      if (hasAcceptableSize) {
        QmmmModel model = {qmRegion, indices, molecularCharge, spinMultiplicity};
        qmmmModelCandidates.push_back(model);
        allQmAtomIndices.push_back(indices);
      }
      if (isReferenceModel) {
        QmmmModel model = {qmRegion, indices, molecularCharge, spinMultiplicity};
        qmmmReferenceModels.push_back(model);
        allQmAtomIndices.push_back(indices);
      }
    }
    r += 0.1;
  }
  // Get rid of ref models if there are too many
  int numRefModelsToErase = qmmmReferenceModels.size() - settings.getInt(maxNumRefModels);
  if (numRefModelsToErase > 0)
    qmmmReferenceModels.erase(qmmmReferenceModels.begin(), qmmmReferenceModels.begin() + numRefModelsToErase);
}

void QmRegionCandidateGenerator::generateSingleQmRegionNonRandomly(std::vector<QmmmModel>& qmmmModelCandidates,
                                                                   SwooseUtilities::SubsystemGenerator& generator,
                                                                   const SwooseUtilities::FragmentAnalyzer& fragmentAnalyzer,
                                                                   const double& initialRadius, int centerAtom) {
  std::vector<int> indices;
  Utils::AtomCollection qmRegion =
      generator.generateSubsystem(centerAtom, indices, initialRadius * Utils::Constants::bohr_per_angstrom);
  std::sort(indices.begin(), indices.end());
  QmmmModel model = {qmRegion, indices, fragmentAnalyzer.getMolecularCharge(), fragmentAnalyzer.getSpinMultiplicity()};
  qmmmModelCandidates.push_back(model);
}

} // namespace Qmmm
} // namespace Scine