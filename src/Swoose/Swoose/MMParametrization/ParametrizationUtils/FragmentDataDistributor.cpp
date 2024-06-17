/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "FragmentDataDistributor.h"
#include "../ParametrizationData.h"
#include <Core/Log.h>

namespace Scine {
namespace MMParametrization {

FragmentDataDistributor::FragmentDataDistributor(ParametrizationData& data, Core::Log& log) : data_(data), log_(log) {
}

bool FragmentDataDistributor::referenceDataIsSufficient(bool refineConnectivity, std::vector<int> failedCalculations,
                                                        std::vector<std::vector<std::string>> failedCalculationIds) const {
  // First, handle the case of no fragmentation
  if (data_.vectorOfStructures.size() == 1) {
    std::vector<int> candidates = {0};
    return fragmentIsCoveredByData(0, candidates, refineConnectivity, failedCalculations, failedCalculationIds);
  }

  // This is for the case when fragmentation took place
  bool uncoveredFragmentFound = false;

  for (int i = 0; i < int(data_.vectorOfStructures.size()); ++i) {
    // Ignore, if superfluous.
    if (std::find(data_.superfluousFragments.begin(), data_.superfluousFragments.end(), i) != data_.superfluousFragments.end())
      continue;

    std::vector<int> candidates = getCandidateFragments(i);
    if (candidates.size() < 4) // important for terminal atoms
      updateCandidateFragmentsWithThirdShellNeighbors(i, candidates);
    if (!fragmentIsCoveredByData(i, candidates, refineConnectivity, failedCalculations, failedCalculationIds))
      uncoveredFragmentFound = true;
  }
  return !uncoveredFragmentFound;
}

bool FragmentDataDistributor::fragmentIsCoveredByData(int fragmentIndex, const std::vector<int>& candidates,
                                                      bool refineConnectivity, const std::vector<int>& failedCalculations,
                                                      std::vector<std::vector<std::string>> failedCalculationIds) const {
  // Check whether the candidate is covered already by the available reference data
  for (const auto& candidate : candidates) {
    if (data_.vectorOfHessians.at(candidate) == nullptr) // If Hessian worked, then the structure opt. also worked.
      continue;
    if (data_.atomicChargesForEachFragment.at(candidate).empty())
      continue;
    if (refineConnectivity && data_.vectorOfBondOrderCollections.at(candidate) == nullptr)
      continue;
    return true;
  }

  // Check whether this fragment is hopeless. If it is, throw exception.
  if (fragmentIsHopeless(candidates, failedCalculations)) {
    log_.debug << "For fragment " << fragmentIndex
               << " all reference calculations failed. Please check the following calculations :" << Core::Log::endl;
    auto candidateCalculations = failedCalculationIds.at(fragmentIndex);
    for (auto candidate : candidateCalculations)
      log_.debug << candidate << Core::Log::endl;

    throw std::runtime_error("Too many calculations failed. Parametrization cannot be successful anymore.");
  }
  return false;
}

bool FragmentDataDistributor::fragmentIsHopeless(const std::vector<int>& candidates,
                                                 const std::vector<int>& failedCalculations) const {
  if (failedCalculations.size() != data_.vectorOfStructures.size())
    return false;

  std::array<int, 3> failsForEach = {0, 0, 0};
  for (const auto& candidate : candidates) {
    int score = failedCalculations.at(candidate);
    if (score % 2 == 0)
      failsForEach[0]++;
    if (score % 3 == 0)
      failsForEach[1]++;
    if (score % 5 == 0)
      failsForEach[2]++;
  }

  for (int i = 0; i < 3; ++i) {
    if (failsForEach[i] == int(candidates.size()))
      return true;
  }
  return false;
}

std::vector<int> FragmentDataDistributor::getCandidateFragments(int fragmentIndex) const {
  std::vector<int> listOfCandidates = {};
  updateCandidateFragments(fragmentIndex, listOfCandidates);
  return listOfCandidates;
}

void FragmentDataDistributor::updateCandidateFragments(int fragmentIndex, std::vector<int>& listOfCandidates) const {
  // Trivial first candidate:
  if (std::find(listOfCandidates.begin(), listOfCandidates.end(), fragmentIndex) == listOfCandidates.end())
    listOfCandidates.push_back(fragmentIndex);
  // For loop over all neighbors and neighbors of neighbors:
  for (const auto& neighbor : data_.listsOfNeighbors.at(fragmentIndex)) {
    if (std::find(listOfCandidates.begin(), listOfCandidates.end(), neighbor) == listOfCandidates.end())
      listOfCandidates.push_back(neighbor);
    for (const auto& secondNeighbor : data_.listsOfNeighbors.at(neighbor)) {
      if (std::find(listOfCandidates.begin(), listOfCandidates.end(), secondNeighbor) == listOfCandidates.end())
        listOfCandidates.push_back(secondNeighbor);
    }
  }
}

void FragmentDataDistributor::updateCandidateFragmentsWithThirdShellNeighbors(int fragmentIndex,
                                                                              std::vector<int>& listOfCandidates) const {
  for (const auto& n1 : data_.listsOfNeighbors.at(fragmentIndex)) {
    for (const auto& n2 : data_.listsOfNeighbors.at(n1)) {
      for (const auto& n3 : data_.listsOfNeighbors.at(n2)) {
        if (std::find(listOfCandidates.begin(), listOfCandidates.end(), n3) == listOfCandidates.end())
          listOfCandidates.push_back(n3);
      }
    }
  }
}

} // namespace MMParametrization
} // namespace Scine
