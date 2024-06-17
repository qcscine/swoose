/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmmmModelAnalyzer.h"
#include "QmRegionSelector.h"
#include "QmRegionSelectorSettings.h"
#include "QmmmReferenceDataManager.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Settings.h>

namespace Scine {
namespace Qmmm {

QmmmModelAnalyzer::QmmmModelAnalyzer(const Utils::Settings& settings, Core::Log& log, const QmmmData& data,
                                     const Utils::AtomCollection& structure, const std::vector<QmmmModel>& qmmmModelCandidates)
  : settings_(settings), log_(log), data_(data), structure_(structure), candidates_(qmmmModelCandidates) {
  if (candidates_.empty())
    throw std::runtime_error("The list of QM/MM model candidates that was given to the analyzer is empty.");
  analyzeData();
}

int QmmmModelAnalyzer::getIndexOfOptimalModel() const {
  return indexOfOptimalModel_;
}

void QmmmModelAnalyzer::analyzeData() {
  const double tolErrorPercentage = settings_.getDouble(SwooseUtilities::SettingsNames::tolerancePercentageError);
  auto relevantAtoms = getAtomIndicesCloseToCenterAtom();
  auto referenceForces = calculateReferenceForces(relevantAtoms);

  std::vector<double> errors;
  std::vector<int> numLinkAtoms;
  for (int i = 0; i < int(candidates_.size()); ++i) {
    auto error = calculateMeanErrorForCandidateModel(i, relevantAtoms, referenceForces);
    errors.push_back(error);
    numLinkAtoms.push_back(data_.linkAtomNumbers.at(i));
  }
  double minimumError = *std::min_element(errors.begin(), errors.end());
  if (minimumError > (std::numeric_limits<double>::max() * 1e-1))
    throw std::runtime_error("None of the calculations for the candidate models were completed successfully.");
  double toleranceError = std::max((tolErrorPercentage / 100.0) * minimumError, 1e-4);
  int minimumNumLinkAtoms = *std::min_element(numLinkAtoms.begin(), numLinkAtoms.end());

  int currentOptimalModelIndex = -1;
  while (currentOptimalModelIndex == -1) {
    for (int j = 0; j < int(errors.size()); ++j) {
      if (numLinkAtoms.at(j) != minimumNumLinkAtoms) {
        continue;
      }
      if (errors.at(j) > (minimumError + toleranceError)) {
        continue;
      }
      if (currentOptimalModelIndex == -1) {
        currentOptimalModelIndex = j;
        continue;
      }
      if (std::max(errors.at(j), toleranceError) < std::max(errors.at(currentOptimalModelIndex), toleranceError))
        currentOptimalModelIndex = j;
    }
    // If no models are within the error threshold with the minimum number of link atoms, increase this number:
    minimumNumLinkAtoms++;
  }
  indexOfOptimalModel_ = currentOptimalModelIndex;
}

std::vector<int> QmmmModelAnalyzer::getAtomIndicesCloseToCenterAtom() {
  auto centerAtoms = settings_.getIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms);
  std::vector<int> relevantAtoms;
  std::vector<Eigen::RowVector3d> distVecs;
  std::vector<int> checkedElements;
  for (int i = 0; i < structure_.size(); ++i) {
    for (auto& centerAtom : centerAtoms) {
      Eigen::RowVector3d distVec = structure_.getPosition(i) - structure_.getPosition(centerAtom);
      // TODO: Reconsider also including hydrogens here
      if (distVec.norm() <= distanceThresholdForAnalysis_ && (Utils::ElementInfo::Z(structure_.getElement(i)) > 1)) {
        // Only add atom once
        if (std::find(relevantAtoms.begin(), relevantAtoms.end(), i) == relevantAtoms.end())
          relevantAtoms.push_back(i);
      }
    }
  }
  return relevantAtoms;
}

std::vector<Eigen::RowVector3d> QmmmModelAnalyzer::calculateReferenceForces(const std::vector<int>& relevantAtoms) {
  std::vector<Eigen::RowVector3d> referenceForces;
  referenceForces.reserve(relevantAtoms.size());
  int refStartIndex = data_.forces.size() - data_.nRef;

  // Loop over all relevant atoms
  for (const auto& a : relevantAtoms) {
    Eigen::RowVector3d force;
    force.setZero();
    int numFailures = 0;
    // Loop over reference systems
    for (int i = refStartIndex; i < int(data_.forces.size()); ++i) {
      if (data_.forces.at(i).rows() != structure_.size()) {
        numFailures++;
        continue;
      }
      Eigen::RowVector3d f = data_.forces.at(i).row(a);
      force += f;
    }
    if (numFailures == data_.nRef)
      throw std::runtime_error("None of the calculations for the reference models were completed successfully.");
    force /= (data_.nRef - numFailures);
    referenceForces.push_back(force);
  }
  return referenceForces;
}

double QmmmModelAnalyzer::calculateMeanErrorForCandidateModel(int modelIndex, const std::vector<int>& relevantAtoms,
                                                              const std::vector<Eigen::RowVector3d>& referenceForces) {
  assert(candidates_.size() == (data_.forces.size() - data_.nRef));
  assert(relevantAtoms.size() == referenceForces.size());
  ForcesCollection forcesForCandidate = data_.forces.at(modelIndex);
  if (forcesForCandidate.rows() != structure_.size())
    return std::numeric_limits<double>::max();
  double totalError = 0.0;

  // Loop over all relevant atoms
  for (int i = 0; i < int(relevantAtoms.size()); ++i) {
    Eigen::RowVector3d diff = forcesForCandidate.row(relevantAtoms[i]) - referenceForces.at(i);
    double error = (std::abs(diff.x()) + std::abs(diff.y()) + std::abs(diff.z())) / 3;
    totalError += error;
  }
  return totalError / relevantAtoms.size();
}

} // namespace Qmmm
} // namespace Scine
