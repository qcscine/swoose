/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AtomicChargesAssembler.h"
#include "../ParametrizationData.h"
#include "FragmentDataDistributor.h"
#include <numeric>

namespace Scine {
namespace MMParametrization {
namespace AtomicChargesAssembler {

void assembleAtomicCharges(ParametrizationData& data, Core::Log& log) {
  int numberOfStructures = data.vectorOfStructures.size();
  data.atomicCharges.resize(data.numberOfAtoms);

  if (int(data.atomicChargesForEachFragment.size()) != numberOfStructures) {
    std::runtime_error(
        "Failure while assembling the atomic charges vector for the whole system. The reference data is incomplete.");
  }

  if (numberOfStructures == 1 && data.atomicChargesForEachFragment[0].empty()) {
    throw std::runtime_error("The atomic charges could not be obtained.");
  }
  else if (numberOfStructures == 1) {
    data.atomicCharges = data.atomicChargesForEachFragment[0];
    return;
  }

  // Get candidate fragments from which to get the data from
  FragmentDataDistributor fragmentDataDistributor(data, log);
  for (int i = 0; i < numberOfStructures; ++i) {
    // Generate a list of candidate fragments to obtain the atomic charge of atom i from.
    std::vector<int> listOfCandidates = fragmentDataDistributor.getCandidateFragments(i);
    if (listOfCandidates.size() < 4) // important for terminal atoms
      fragmentDataDistributor.updateCandidateFragmentsWithThirdShellNeighbors(i, listOfCandidates);

    // Check whether the calculations were successful for any of the candidate fragments
    bool successful = false;
    for (const auto& candidate : listOfCandidates) {
      if (!data.atomicChargesForEachFragment[candidate].empty()) {
        auto atomIndexInSuccessfullyCalculatedFragment =
            std::distance(data.atomIndexMapping[candidate].begin(),
                          std::find(data.atomIndexMapping[candidate].begin(), data.atomIndexMapping[candidate].end(), i));

        // Throw exception in the very unlikely case of the atom not being present in the candidate fragment
        if (atomIndexInSuccessfullyCalculatedFragment >= int(data.atomIndexMapping[candidate].size()))
          throw std::runtime_error("Error while assembling the atomic charges: The atom of interest is not present "
                                   "in the candidate fragment.");

        // Fill the atomic charges in the ParametrizationData object accordingly
        data.atomicCharges[i] = data.atomicChargesForEachFragment[candidate][atomIndexInSuccessfullyCalculatedFragment];
        successful = true;
        break; // No need to look at the other neighbors in a successful case
      }
    }
    if (!successful) {
      // The following function tries to update this charge with user-provided parameters and returns whether
      // such parameters were available. If not, we throw an exception.
      if (!updateChargeWithExistingParameters(data.atomicCharges[i], i, data))
        throw std::runtime_error("The atomic charge of atom " + std::to_string(i) + " could not be obtained.");
    }
  }
}

void renormalizeAtomicCharges(ParametrizationData& data) {
  // Get sum of atomic charges
  double sumOfAtomicCharges = std::accumulate(data.atomicCharges.begin(), data.atomicCharges.end(), 0.0);

  // Get total charge of the full system from ParametrizationData object
  int totalCharge = 0;
  for (const auto& entry : data.formalCharges)
    totalCharge += entry.second;

  // Renormalize the charges
  double correction = (static_cast<double>(totalCharge) - sumOfAtomicCharges) / data.numberOfAtoms;
  for (auto& charge : data.atomicCharges)
    charge += correction;
}

bool updateChargeWithExistingParameters(double& atomicCharge, int atomIndex, const ParametrizationData& data) {
  std::string atomType = data.atomTypes.getAtomType(atomIndex);
  if (data.parameters.getCharges().find(atomType) != data.parameters.getCharges().end()) {
    atomicCharge = data.parameters.getCharges().at(atomType);
    return true;
  }
  return false;
}

} // namespace AtomicChargesAssembler
} // namespace MMParametrization
} // namespace Scine
