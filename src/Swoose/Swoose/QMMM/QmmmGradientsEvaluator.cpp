/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmmmGradientsEvaluator.h"
#include <Utils/Geometry/AtomCollection.h>

namespace Scine {
namespace Qmmm {

QmmmGradientsEvaluator::QmmmGradientsEvaluator(const Utils::GradientCollection& qmGradients,
                                               const Utils::GradientCollection& mmGradients,
                                               const Utils::GradientCollection& pcGradients,
                                               const std::vector<int>& listOfQmAtoms, const std::vector<int>& mmBoundaryAtoms,
                                               const std::vector<std::list<int>>& listsOfNeighbors,
                                               const Utils::AtomCollection& fullStructure, const Utils::AtomCollection& qmRegion)
  : qmGradients_(qmGradients),
    mmGradients_(mmGradients),
    pcGradients_(pcGradients),
    listOfQmAtoms_(listOfQmAtoms),
    mmBoundaryAtoms_(mmBoundaryAtoms),
    listsOfNeighbors_(listsOfNeighbors),
    fullStructure_(fullStructure),
    qmRegion_(qmRegion) {
}

Utils::GradientCollection QmmmGradientsEvaluator::calculateQmmmGradients() {
  Utils::GradientCollection combinedGradients = mmGradients_;
  for (int i = 0; i < listOfQmAtoms_.size(); ++i) {
    // This works since the link atoms in the QM region are ALWAYS listed at the end of the AtomCollection.
    combinedGradients.row(listOfQmAtoms_.at(i)) += qmGradients_.row(i);
  }

  if (pcGradients_.size() > 0) {
    int pcGradientsRow = 0;
    for (int i = 0; i < combinedGradients.rows(); ++i) {
      if (std::find(listOfQmAtoms_.begin(), listOfQmAtoms_.end(), i) == listOfQmAtoms_.end()) {
        combinedGradients.row(i) += pcGradients_.row(pcGradientsRow);
        pcGradientsRow++;
      }
    }
  }

  addBoundaryGradientsContributions(combinedGradients);

  return combinedGradients;
}

void QmmmGradientsEvaluator::addBoundaryGradientsContributions(Utils::GradientCollection& qmmmGradients) {
  auto linkAtomIndex = listOfQmAtoms_.size(); // First link atom is right after the real QM atoms.
  int mmAtomCounter = 0;

  for (int i = 0; i < listOfQmAtoms_.size(); ++i) {
    int qmAtomIndex = listOfQmAtoms_.at(i);
    auto neighbors = listsOfNeighbors_.at(qmAtomIndex);
    // Iterate over all atoms bonded to the QM atom
    for (const auto& neighbor : neighbors) {
      // Check whether this neighbor is NOT a QM atom
      if (std::find(listOfQmAtoms_.begin(), listOfQmAtoms_.end(), neighbor) == listOfQmAtoms_.end()) {
        // Get the link atom position and increment the link atom index counter.
        Eigen::RowVector3d linkAtomPosition = qmRegion_.getPosition(linkAtomIndex);
        Eigen::RowVector3d linkAtomGradient = qmGradients_.row(linkAtomIndex);
        linkAtomIndex++;

        // Get the position of the corresponding MM atom and increment the MM atom index counter.
        int mmAtomIndex = mmBoundaryAtoms_.at(mmAtomCounter);
        Eigen::RowVector3d mmAtomPosition = fullStructure_.getPosition(mmAtomIndex);
        mmAtomCounter++;

        // Get the position of the corresponding QM atom
        Eigen::RowVector3d qmAtomPosition = fullStructure_.getPosition(qmAtomIndex);

        // Add contributions for this link atom
        for (int k = 0; k < 3; ++k) {
          auto contributions = calculateGradientContributionForOneBoundary(qmAtomPosition, mmAtomPosition,
                                                                           linkAtomPosition, linkAtomGradient, k);
          qmmmGradients(qmAtomIndex, k) += contributions.first;
          qmmmGradients(mmAtomIndex, k) += contributions.second;
        }
      }
    }
  }
}

std::pair<double, double> QmmmGradientsEvaluator::calculateGradientContributionForOneBoundary(
    const Eigen::Ref<Eigen::RowVector3d> qmAtomPosition, const Eigen::Ref<Eigen::RowVector3d> mmAtomPosition,
    const Eigen::Ref<Eigen::RowVector3d> linkAtomPosition, const Eigen::Ref<Eigen::RowVector3d> linkAtomGradient,
    int dimension) {
  Eigen::RowVector3d distanceVec = mmAtomPosition - qmAtomPosition;
  double oneDimDist = mmAtomPosition(dimension) - qmAtomPosition(dimension);

  Eigen::RowVector3d unitVec(0.0, 0.0, 0.0);
  unitVec(dimension) = 1.0;

  double d1 = distanceVec.norm();
  double d2 = (qmAtomPosition - linkAtomPosition).norm();

  Eigen::RowVector3d qmCoordDeriv = (1 - d2 / d1) * unitVec + ((d2 * oneDimDist) / std::pow(d1, 3)) * distanceVec;
  double qmContribution = linkAtomGradient * qmCoordDeriv.transpose();

  Eigen::RowVector3d mmCoordDeriv = (d2 / d1) * unitVec - ((d2 * oneDimDist) / std::pow(d1, 3)) * distanceVec;
  double mmContribution = linkAtomGradient * mmCoordDeriv.transpose();

  return std::make_pair(qmContribution, mmContribution);
}

} // namespace Qmmm
} // namespace Scine