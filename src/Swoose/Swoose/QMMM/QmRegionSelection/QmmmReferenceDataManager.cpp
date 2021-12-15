/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmmmReferenceDataManager.h"
#include "QmRegionSelector.h"
#include "QmRegionSelectorSettings.h"
#include "QmmmDirectCalculationsHelper.h"
#include "SymmetryScores.h"
#include <Core/Log.h>
#include <Utils/Bonds/BondOrderCollection.h>

#ifdef SWOOSE_COMPILE_DATABASE
#  include "QmmmDatabaseHelper.h"
#endif

namespace Scine {
namespace Qmmm {

QmmmReferenceDataManager::QmmmReferenceDataManager(const Utils::Settings& settings, Core::Log& log,
                                                   const Utils::AtomCollection& structure,
                                                   const Utils::BondOrderCollection& bondOrders,
                                                   const std::vector<QmmmModel>& qmmmModelCandidates,
                                                   const std::vector<QmmmModel>& qmmmReferenceModels)
  : settings_(settings),
    log_(log),
    structure_(structure),
    qmmmModelCandidates_(qmmmModelCandidates),
    qmmmReferenceModels_(qmmmReferenceModels),
    bondOrders_(bondOrders) {
}

QmmmData QmmmReferenceDataManager::calculateData() {
  QmmmData data;
  data.nRef = qmmmReferenceModels_.size();
  calculateLinkAtomNumbers(data);
  calculateSymmetryScores(data);
  handleReferenceCalculations(data);
  return data;
}

void QmmmReferenceDataManager::calculateLinkAtomNumbers(QmmmData& data) const {
  data.linkAtomNumbers.clear();
  data.linkAtomNumbers.reserve(qmmmModelCandidates_.size() + qmmmReferenceModels_.size());
  for (const auto& candidateModel : qmmmModelCandidates_)
    data.linkAtomNumbers.push_back(std::count(candidateModel.qmAtomIndices.begin(), candidateModel.qmAtomIndices.end(), -1));
  for (const auto& refModel : qmmmReferenceModels_)
    data.linkAtomNumbers.push_back(std::count(refModel.qmAtomIndices.begin(), refModel.qmAtomIndices.end(), -1));
}

void QmmmReferenceDataManager::calculateSymmetryScores(QmmmData& data) const {
  data.symmetryScores.clear();
  data.symmetryScores.reserve(qmmmModelCandidates_.size() + qmmmReferenceModels_.size());
  int centerAtom = settings_.getInt(SwooseUtilities::SettingsNames::qmRegionCenterAtom);

  // Pre-calculate the distances of each atom to the center atom only once:
  std::vector<double> distances(structure_.size());
  for (int i = 0; i < structure_.size(); ++i) {
    Eigen::RowVector3d distVec = structure_.getPosition(i) - structure_.getPosition(centerAtom);
    distances[i] = distVec.norm();
  }

  // Calculate symmetry scores for candidate models
  for (const auto& candidateModel : qmmmModelCandidates_) {
    auto score = calculateSymmetryScoreForOneModel(candidateModel, distances, 3);
    data.symmetryScores.push_back(score);
  }

  // Calculate symmetry scores for reference models
  for (const auto& refModel : qmmmReferenceModels_) {
    auto score = calculateSymmetryScoreForOneModel(refModel, distances, 3);
    data.symmetryScores.push_back(score);
  }
}

void QmmmReferenceDataManager::handleReferenceCalculations(QmmmData& data) {
  auto mode = settings_.getString(SwooseUtilities::SettingsNames::referenceDataMode);

  if (mode == SwooseUtilities::OptionNames::databaseMode) {
#ifdef SWOOSE_COMPILE_DATABASE
    QmmmDatabaseHelper databaseHelper(settings_, log_, structure_, bondOrders_, qmmmModelCandidates_,
                                      qmmmReferenceModels_, data);
    data.forces = databaseHelper.calculateForces();
#else
    throw std::runtime_error("Swoose was not compiled with database support.");
#endif
  }
  else if (mode == SwooseUtilities::OptionNames::directMode) {
    QmmmDirectCalculationsHelper directCalculationsHelper(settings_, log_, structure_, qmmmModelCandidates_,
                                                          qmmmReferenceModels_, data);
    data.forces = directCalculationsHelper.calculateForces();
  }
  else {
    throw std::runtime_error(
        "Currently, only database mode and direct mode are implemented for calculating QM/MM model forces.");
  }
}

} // namespace Qmmm
} // namespace Scine