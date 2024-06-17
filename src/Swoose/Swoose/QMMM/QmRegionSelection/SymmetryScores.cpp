/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SymmetryScores.h"
#include "QmRegionSelector.h"
#include "QmRegionSelectorSettings.h"
#include "QmmmReferenceDataManager.h"

namespace Scine {
namespace Qmmm {

double calculateSymmetryScoreForOneModel(const QmmmModel& model, const std::vector<double>& distances, int numberOfValuesForMean) {
  std::vector<double> distancesOfIncludedAtoms;
  std::vector<double> distancesOfExcludedAtoms;
  for (int i = 0; i < int(distances.size()); ++i) {
    if (std::find(model.qmAtomIndices.begin(), model.qmAtomIndices.end(), i) != model.qmAtomIndices.end())
      distancesOfIncludedAtoms.push_back(distances[i]);
    else
      distancesOfExcludedAtoms.push_back(distances[i]);
  }

  std::sort(distancesOfExcludedAtoms.begin(), distancesOfExcludedAtoms.end());
  std::sort(distancesOfIncludedAtoms.begin(), distancesOfIncludedAtoms.end(), std::greater<>());

  double meanOfSmallest = 0.0;
  double meanOfLargest = 0.0;
  try {
    for (int k = 0; k < numberOfValuesForMean; ++k) {
      meanOfSmallest += distancesOfExcludedAtoms.at(k);
    }
    meanOfSmallest /= numberOfValuesForMean;
    for (int l = 0; l < numberOfValuesForMean; ++l) {
      meanOfLargest += distancesOfIncludedAtoms.at(l);
    }
    meanOfLargest /= numberOfValuesForMean;
  }
  catch (const std::out_of_range& e) {
    return 0.0;
  }
  return meanOfLargest / meanOfSmallest;
}

double calculateMaxAllowedSymmetryScore(const QmmmData& qmmmData, const Utils::Settings& settings) {
  if (qmmmData.symmetryScores.empty())
    throw std::runtime_error("Symmetry scores could not be evaluated.");

  if (settings.getIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms).size() > 1) {
    std::cout << "If you choose multiple QM center atoms, the symmetry of the generated overall QM region will not "
                 "be considered."
              << std::endl;
  }

  // Evaluate symmetry score only for single-atom QM regions
  if (settings.getIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms).size() == 1) {
    const double minSymmetryScore =
        *std::min_element(qmmmData.symmetryScores.begin(), qmmmData.symmetryScores.end() - qmmmData.nRef);
    const double tolerance =
        (settings.getDouble(SwooseUtilities::SettingsNames::tolerancePercentageSymmetryScore) / 100.0) * minSymmetryScore;
    return minSymmetryScore + tolerance;
  }
  else {
    const double maxSymmetryScore =
        *std::max_element(qmmmData.symmetryScores.begin(), qmmmData.symmetryScores.end() - qmmmData.nRef);
    return maxSymmetryScore;
  }
}

} // namespace Qmmm
} // namespace Scine
