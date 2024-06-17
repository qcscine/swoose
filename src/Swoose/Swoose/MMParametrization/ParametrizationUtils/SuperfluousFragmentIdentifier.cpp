/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SuperfluousFragmentIdentifier.h"
#include "../ParametrizationData.h"
#include "ReparametrizationHelper.h"
#include <Core/Log.h>

namespace Scine {
namespace MMParametrization {
namespace SuperfluousFragmentIdentifier {

void identifySuperfluousFragments(ParametrizationData& data, Core::Log& log,
                                  std::shared_ptr<ReparametrizationHelper> reparametrizationHelper) {
  data.superfluousFragments.clear(); // should already be empty, but for consistency.

  if (data.vectorOfStructures.size() == 1)
    return;

  for (int i = 0; i < int(data.vectorOfStructures.size()); ++i) {
    // It may be superfluous since we are in a process of local re-parametrization of the model:
    if (reparametrizationHelper != nullptr) { // Important! It is null if we are not in such a process.
      if (!reparametrizationHelper->isRelevantFragment(i) && !data.siteIspHSensitive.at(i)) {
        data.superfluousFragments.push_back(i);
        continue;
      }
    }

    // if center atom is pH sensitive, this structure should not be superfluous
    if (data.pHSensitiveSites.count(i) == 1)
      continue;

    // Check for duplicate structures:
    Utils::AtomCollection structure = *data.vectorOfStructures.at(i);
    auto charge = data.vectorOfChargesAndMultiplicities.at(i).first;
    auto mult = data.vectorOfChargesAndMultiplicities.at(i).second;
    auto indices = data.atomIndexMapping.at(i);

    bool original = true;
    int correspondingOriginal = -1;
    for (int j = 0; j < i; ++j) {
      Utils::AtomCollection otherStructure = *data.vectorOfStructures.at(j);

      // Check 1: number of atoms
      if (otherStructure.size() != structure.size())
        continue;

      // Check 2: molecular charge
      if (charge != data.vectorOfChargesAndMultiplicities.at(j).first)
        continue;

      // Check 3: spin multiplicity
      if (mult != data.vectorOfChargesAndMultiplicities.at(j).second)
        continue;

      // Check 4: all atoms are the same
      bool allAtomsInside = true;
      auto otherIndices = data.atomIndexMapping.at(j);
      for (int index : indices) {
        if (std::find(otherIndices.begin(), otherIndices.end(), index) == otherIndices.end()) {
          allAtomsInside = false;
          break;
        }
      }

      // Decision:
      if (allAtomsInside) {
        if (std::find(data.superfluousFragments.begin(), data.superfluousFragments.end(), j) ==
            data.superfluousFragments.end()) { // original is only original if it wasn't removed yet
          original = false;
          correspondingOriginal = j;
          break;
        }
      }
    }

    // If structure is not an original
    if (!original) {
      // Check whether the corresponding original is a neighboring atom
      if (std::find(data.listsOfNeighbors.at(i).begin(), data.listsOfNeighbors.at(i).end(), correspondingOriginal) !=
          data.listsOfNeighbors.at(i).end()) {
        // Add this fragment index to the data object
        data.superfluousFragments.push_back(i);
      }
    }
  }

  log.output << "Number of superfluous fragments: " << data.superfluousFragments.size() << Core::Log::endl;
  log.debug << "Superfluous fragments: ";
  for (const auto& f : data.superfluousFragments)
    log.debug << f << " ";
  log.debug << Core::Log::endl;
}

} // namespace SuperfluousFragmentIdentifier
} // namespace MMParametrization
} // namespace Scine
