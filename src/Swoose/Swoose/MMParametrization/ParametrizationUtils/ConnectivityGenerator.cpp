/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ConnectivityGenerator.h"
#include "../MMParametrizationSettings.h"
#include "../ParametrizationData.h"
#include "FragmentDataDistributor.h"
#include <Core/Log.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondDetector.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <boost/filesystem.hpp>
#include <chrono>

namespace Scine {
namespace MMParametrization {

ConnectivityGenerator::ConnectivityGenerator(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings,
                                             Core::Log& log)
  : data_(data), settings_(settings), log_(log) {
  bondOrderThreshold_ = settings_->getDouble(SwooseUtilities::SettingsNames::bondOrderThreshold);
}

void ConnectivityGenerator::generateInitialListsOfNeighbors() {
  // If connectivity file exists, use that
  auto existingConnFile = settings_->getString(SwooseUtilities::SettingsNames::connectivityFilePath);
  if (!existingConnFile.empty() && boost::filesystem::exists(existingConnFile)) {
    log_.output << "Reading connectivity from file: " << existingConnFile << Core::Log::endl;
    data_.listsOfNeighbors = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(existingConnFile);

    // Check if connectivity file was valid
    if (data_.listsOfNeighbors.size() != data_.numberOfAtoms)
      throw std::runtime_error(
          "The number of atoms in the provided connectivity file does not match the one of the molecular structure.");

    data_.bondOrders = SwooseUtilities::TopologyUtils::generateBondOrderMatrixFromListsOfNeighbors(data_.listsOfNeighbors);
    return;
  }
  // Bond detection with covalent radii
  data_.bondOrders = Utils::BondDetector::detectBonds(data_.fullStructure);
  data_.listsOfNeighbors =
      SwooseUtilities::TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(data_.numberOfAtoms, data_.bondOrders, 0.5);
  log_.output << "Obtained the system connectivity from covalent radii bond detector." << Core::Log::endl;
}

void ConnectivityGenerator::refineListsOfNeighbors() {
  assert(data_.vectorOfBondOrderCollections.size() == data_.vectorOfStructures.size());
  std::vector<std::list<int>> initialListsOfNeighbors = data_.listsOfNeighbors;

  // Initialize for later use
  FragmentDataDistributor fragmentDataDistributor(data_);
  // Loop over all atoms
  std::vector<std::vector<int>> candidates;
  candidates.reserve(data_.numberOfAtoms);
  for (int i = 0; i < data_.numberOfAtoms; ++i) {
    // Get candidate atoms for using their bond order matrix
    std::vector<int> candidatesForBondOrderMatrix = {0}; // stays like this if the statement below is false
    if (data_.vectorOfStructures.size() > 1) {
      candidatesForBondOrderMatrix = fragmentDataDistributor.getCandidateFragments(i);
      if (candidatesForBondOrderMatrix.size() < 4) // important for terminal atoms
        fragmentDataDistributor.updateCandidateFragmentsWithThirdShellNeighbors(i, candidatesForBondOrderMatrix);
    }
    candidates.push_back(candidatesForBondOrderMatrix);
  }

  // Clear the previous bond orders
  for (int i = 0; i < data_.numberOfAtoms; ++i) {
    data_.listsOfNeighbors[i].clear();
  }

  int warningCounter = 0;
  for (int i = 0; i < data_.numberOfAtoms; ++i) {
    bool bondOrdersObtained = false;
    // Loop over all candidates until a bond order matrix is available
    for (const auto& candidate : candidates.at(i)) {
      // Check whether the bond order matrix for that candidate is available
      if (data_.vectorOfBondOrderCollections[candidate]) {
        bondOrdersObtained = true;
        // Loop over all atoms
        for (int j = 0; j < data_.numberOfAtoms; ++j) {
          // Get indices of the atoms i and j in the candidate fragment
          auto indexOfI =
              std::distance(data_.atomIndexMapping[candidate].begin(),
                            std::find(data_.atomIndexMapping[candidate].begin(), data_.atomIndexMapping[candidate].end(), i));
          auto indexOfJ =
              std::distance(data_.atomIndexMapping[candidate].begin(),
                            std::find(data_.atomIndexMapping[candidate].begin(), data_.atomIndexMapping[candidate].end(), j));

          // Make sure both atoms are actually in the given fragment
          auto numberOfAtomsInFragments = data_.atomIndexMapping[candidate].size();
          if ((indexOfI >= numberOfAtomsInFragments) || (indexOfJ >= numberOfAtomsInFragments))
            continue;

          // Check the bond order of i and j
          if (data_.vectorOfBondOrderCollections[candidate]->getOrder(indexOfI, indexOfJ) > bondOrderThreshold_) {
            if (std::find(data_.listsOfNeighbors[i].begin(), data_.listsOfNeighbors[i].end(), j) ==
                data_.listsOfNeighbors[i].end()) {
              data_.listsOfNeighbors[i].push_back(j);
            }
            if (std::find(data_.listsOfNeighbors[j].begin(), data_.listsOfNeighbors[j].end(), i) ==
                data_.listsOfNeighbors[j].end()) {
              data_.listsOfNeighbors[j].push_back(i);
            }
          }
        }
        // Break after one candidate was successful
        break;
      }
    }
    // If the bond orders could not be obtained, give a warning and use the initial bond orders for compensation
    if (!bondOrdersObtained) {
      warningCounter++;
      for (const auto& neighbor : initialListsOfNeighbors[i]) {
        if (std::find(data_.listsOfNeighbors[i].begin(), data_.listsOfNeighbors[i].end(), neighbor) ==
            data_.listsOfNeighbors[i].end()) {
          data_.listsOfNeighbors[i].push_back(neighbor);
        }
        if (std::find(data_.listsOfNeighbors[neighbor].begin(), data_.listsOfNeighbors[neighbor].end(), i) ==
            data_.listsOfNeighbors[neighbor].end()) {
          data_.listsOfNeighbors[neighbor].push_back(i);
        }
      }
    }
  }
  if (warningCounter > 0) {
    log_.output << "For some atoms, the bond orders could not be obtained. The information from the initial "
                   "connectivity was used for these atoms instead."
                << Core::Log::nl << "Number of atoms with this issue: " << warningCounter << Core::Log::endl;
  }
  // Free memory at the end of this function, we don't need those matrices anymore
  data_.vectorOfBondOrderCollections.clear();
}

} // namespace MMParametrization
} // namespace Scine