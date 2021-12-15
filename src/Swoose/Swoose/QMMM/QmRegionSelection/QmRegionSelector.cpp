/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmRegionSelector.h"
#include "QmRegionCandidateGenerator.h"
#include "QmRegionSelectorSettings.h"
#include "QmmmModelAnalyzer.h"
#include "QmmmReferenceDataManager.h"
#include <Core/Log.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/Geometry/AtomCollection.h>

namespace Scine {
namespace Qmmm {

QmRegionSelector::QmRegionSelector() {
  this->settings_ = std::make_unique<QmRegionSelectorSettings>();
}

void QmRegionSelector::generateQmRegion(const Utils::AtomCollection& fullSystem) {
  // First, some clean-up
  qmmmModelCandidates_.clear();
  qmmmReferenceModels_.clear();
  Utils::BondOrderCollection bondOrders = this->getBondOrders(fullSystem);
  QmRegionCandidateGenerator::generateQmRegionCandidates(qmmmModelCandidates_, qmmmReferenceModels_, fullSystem,
                                                         bondOrders, this->settings(), this->getLog());

  this->getLog().output << "Candidate and reference models were generated." << Core::Log::nl;
  if (qmmmModelCandidates_.size() == 1) {
    this->getLog().output << "Only one QM/MM candidate model was generated." << Core::Log::nl << Core::Log::endl;
    selectedQmRegionIndex_ = 0;
    return;
  }

  this->getLog().output << Core::Log::nl << "Number of QM/MM model candidates: " << qmmmModelCandidates_.size()
                        << Core::Log::nl;
  this->getLog().output << "Number of QM/MM reference models: " << qmmmReferenceModels_.size() << Core::Log::endl
                        << Core::Log::endl;

  if (qmmmModelCandidates_.empty())
    throw std::runtime_error("No QM/MM model candidates could be constructed.");
  if (qmmmReferenceModels_.empty())
    throw std::runtime_error("No QM/MM reference models could be constructed.");

  // Reference data generation if there is more than one candidate
  this->getLog().output << "Starting reference calculations..." << Core::Log::nl << Core::Log::endl;
  QmmmReferenceDataManager referenceDataManager(this->settings(), this->getLog(), fullSystem, bondOrders,
                                                qmmmModelCandidates_, qmmmReferenceModels_);
  auto data = referenceDataManager.calculateData();

  // Analysis
  this->getLog().output << "Analyzing the data to select the optimal QM region..." << Core::Log::nl << Core::Log::endl;
  QmmmModelAnalyzer analyzer(this->settings(), this->getLog(), data, fullSystem, qmmmModelCandidates_);
  selectedQmRegionIndex_ = analyzer.getIndexOfOptimalModel();
}

std::vector<int> QmRegionSelector::getQmRegionIndices() const {
  if (selectedQmRegionIndex_ == -1)
    throw QmRegionHasNotBeenSelectedException();
  return qmmmModelCandidates_.at(selectedQmRegionIndex_).qmAtomIndices;
}

Utils::AtomCollection QmRegionSelector::getQmRegionStructure() const {
  if (selectedQmRegionIndex_ == -1)
    throw QmRegionHasNotBeenSelectedException();
  return qmmmModelCandidates_.at(selectedQmRegionIndex_).structure;
}

std::pair<int, int> QmRegionSelector::getQmRegionChargeAndMultiplicity() const {
  if (selectedQmRegionIndex_ == -1)
    throw QmRegionHasNotBeenSelectedException();
  return {qmmmModelCandidates_.at(selectedQmRegionIndex_).molecularCharge,
          qmmmModelCandidates_.at(selectedQmRegionIndex_).spinMultiplicity};
}

Utils::Settings& QmRegionSelector::settings() {
  return *settings_;
}

const Utils::Settings& QmRegionSelector::settings() const {
  return *settings_;
}

Utils::BondOrderCollection QmRegionSelector::getBondOrders(const Utils::AtomCollection& structure) const {
  auto connFile = settings_->getString(SwooseUtilities::SettingsNames::connectivityFilePath);
  auto listsOfNeighbors = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(connFile);
  if (listsOfNeighbors.size() != structure.size())
    throw std::runtime_error(
        "The number of atoms in the provided connectivity file does not match the one of the molecular structure.");
  return SwooseUtilities::TopologyUtils::generateBondOrderMatrixFromListsOfNeighbors(listsOfNeighbors);
}

} // namespace Qmmm
} // namespace Scine