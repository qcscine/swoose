/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmRegionCandidateGenerator.h"
#include "../QmmmHelpers.h"
#include "QmRegionSelector.h"
#include "QmRegionSelectorSettings.h"
#include <Swoose/Utilities/AtomicInformationReader.h>
#include <Swoose/Utilities/FragmentAnalyzer.h>
#include <Swoose/Utilities/FragmentationHelper.h>
#include <Swoose/Utilities/SubsystemGenerator.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/NativeFilenames.h>
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
  auto centerAtoms = settings.getIntList(qmRegionCenterAtoms);
  int minSize = settings.getInt(qmRegionCandidateMinSize);
  int maxSize = settings.getInt(qmRegionCandidateMaxSize);
  int maxSizeRef = settings.getInt(qmRegionRefMaxSize);
  int minSizeRef = static_cast<int>(0.95 * maxSizeRef); // TODO: Also setting?
  int randomSeed = settings.getInt(qmRegionSelectionRandomSeed);
  auto listsOfNeighbors =
      SwooseUtilities::TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(fullStructure.size(), bondOrders, 0.4);

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
  else {
    log.output
        << "No atomic info file provided. All QM region candidates will be assigned a charge of zero and a spin "
           "multiplicity of one. If you have charged species in your structure, please provide an atomic info file."
        << Core::Log::nl;
  }

  SwooseUtilities::FragmentAnalyzer fragmentAnalyzer(formalCharges, unpairedElectrons);
  SwooseUtilities::SubsystemGenerator generator(fullStructure, bondOrders, fragmentAnalyzer, 0.5,
                                                std::numeric_limits<int>::max(), log, randomSeed, probability);

  if (probability == 1.0) {
    log.output << "Because a cutting probability of 1.0 was selected, only one QM region will be generated based on "
                  "the provided initial radius setting.\nIf the resulting QM region is too large for your purposes, "
                  "re-try with a smaller value for the initial radius."
               << Core::Log::nl << Core::Log::endl;
    auto model = generateSingleQmRegionNonRandomly(generator, fragmentAnalyzer, initialRadius, centerAtoms,
                                                   fullStructure, listsOfNeighbors);
    qmmmModelCandidates.push_back(model);
    return;
  }
  std::vector<std::vector<int>> allQmAtomIndices;
  double r = initialRadius;
  int numExceededRefSizeLimit = 0;
  while (numExceededRefSizeLimit != numStructures) { // Stop if every qm region for given r is too large
    numExceededRefSizeLimit = 0;
    for (int k = 0; k < numStructures; ++k) {
      std::vector<int> indices;
      Utils::AtomCollection qmRegion;
      std::vector<int> mmBoundaryAtoms;
      if (centerAtoms.size() > 1) {
        std::vector<std::vector<int>> listOfMappings;
        std::vector<Utils::AtomCollection> subRegions;
        for (auto& centerAtom : centerAtoms) {
          std::vector<int> tempIndices;
          Utils::AtomCollection tempQmRegion =
              generator.generateSubsystem(centerAtom, tempIndices, initialRadius * Utils::Constants::bohr_per_angstrom);
          listOfMappings.push_back(tempIndices);
          subRegions.push_back(tempQmRegion);
        }
        qmRegion = SwooseUtilities::FragmentationHelper::mergeSubsystems(indices, subRegions, listOfMappings);
        QmmmHelpers::addAllLinkAtoms(qmRegion, fullStructure, listsOfNeighbors, indices, mmBoundaryAtoms);
        indices.insert(indices.end(), mmBoundaryAtoms.size(), -1);
      }
      else
        qmRegion = generator.generateSubsystem(centerAtoms[0], indices, r * Utils::Constants::bohr_per_angstrom);

      std::sort(indices.begin(), indices.end());
      bool valid = fragmentAnalyzer.analyzeFragment(qmRegion, indices);
      if (!valid)
        throw std::runtime_error(
            "The given system with its formal charges and numbers of unpaired electrons is not valid.");

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

QmmmModel QmRegionCandidateGenerator::generateSingleQmRegionNonRandomly(SwooseUtilities::SubsystemGenerator& generator,
                                                                        SwooseUtilities::FragmentAnalyzer& fragmentAnalyzer,
                                                                        const double& initialRadius, std::vector<int> centerAtoms,
                                                                        const Utils::AtomCollection fullStructure,
                                                                        const std::vector<std::list<int>> listsOfNeighbors) {
  std::vector<int> indices;
  Utils::AtomCollection qmRegion;

  if (centerAtoms.size() > 1) {
    std::vector<std::vector<int>> listOfMappings;
    std::vector<Utils::AtomCollection> subRegions;
    std::vector<int> mmBoundaryAtoms;
    for (auto& centerAtom : centerAtoms) {
      std::vector<int> tempIndices;
      Utils::AtomCollection subRegion =
          generator.generateSubsystem(centerAtom, tempIndices, initialRadius * Utils::Constants::bohr_per_angstrom);
      listOfMappings.push_back(tempIndices);
      subRegions.push_back(subRegion);
    }
    qmRegion = SwooseUtilities::FragmentationHelper::mergeSubsystems(indices, subRegions, listOfMappings);
    QmmmHelpers::addAllLinkAtoms(qmRegion, fullStructure, listsOfNeighbors, indices, mmBoundaryAtoms);
    indices.insert(indices.end(), mmBoundaryAtoms.size(), -1);
  }
  else
    qmRegion = generator.generateSubsystem(centerAtoms[0], indices, initialRadius * Utils::Constants::bohr_per_angstrom);

  bool valid = fragmentAnalyzer.analyzeFragment(qmRegion, indices);
  if (!valid)
    throw std::runtime_error(
        "The given system with its formal charges and numbers of unpaired electrons is not valid.");

  std::sort(indices.begin(), indices.end());
  QmmmModel model = {qmRegion, indices, fragmentAnalyzer.getMolecularCharge(), fragmentAnalyzer.getSpinMultiplicity()};
  return model;
}

} // namespace Qmmm
} // namespace Scine
