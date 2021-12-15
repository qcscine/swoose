/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
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
  for (int i = 0; i < distances.size(); ++i) {
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

  const double minSymmetryScore =
      *std::min_element(qmmmData.symmetryScores.begin(), qmmmData.symmetryScores.end() - qmmmData.nRef);
  const double tolerance =
      (settings.getDouble(SwooseUtilities::SettingsNames::tolerancePercentageSymmetryScore) / 100.0) * minSymmetryScore;
  return minSymmetryScore + tolerance;
}

} // namespace Qmmm
} // namespace Scine