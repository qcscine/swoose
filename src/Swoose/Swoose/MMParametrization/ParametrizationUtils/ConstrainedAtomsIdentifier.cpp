/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ConstrainedAtomsIdentifier.h"
#include "../ParametrizationData.h"
#include <Molassembler/Graph.h>
#include <Molassembler/GraphAlgorithms.h>
#include <Molassembler/Interpret.h>
#include <Swoose/Utilities/FragmentationHelper.h>
#include <Utils/Bonds/BondDetector.h>
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace MMParametrization {

ConstrainedAtomsIdentifier::ConstrainedAtomsIdentifier(ParametrizationData& data) : data_(data) {
  data_.constrainedAtoms.resize(data_.numberOfAtoms);
}

void ConstrainedAtomsIdentifier::updateInformationAboutConstrainedAtoms(const Utils::AtomCollection& subsystem,
                                                                        int subsystemIndex) {
  // This will be added to the ParametrizationData object once it is filled
  std::vector<int> constrainedAtoms;

  auto bondOrders = Utils::BondDetector::detectBonds(subsystem);
  auto result = Molassembler::Interpret::graphs(subsystem, bondOrders);
  assert(subsystem.size() == result.componentMap.size());

  // Invert the component map with molassembler
  auto indexMap = result.componentMap.invert();

  // First, gather all of the heavy atoms that are possible candidates for constraining them
  for (int g = 0; g < int(result.graphs.size()); ++g) {
    const auto& graph = result.graphs[g];
    // Get the heavy atoms that were involved in a cleft bond
    std::vector<int> heavyAtomTermini =
        getHeavyAtomTermini(indexMap[g], graph, subsystem, data_.atomIndexMapping[subsystemIndex]);
    // Update the constrainedAtoms vector for this molecule based on the candidates vector (heavy atom termini).
    updateConstrainedAtomsForOneMolecule(constrainedAtoms, heavyAtomTermini, indexMap[g], graph, subsystemIndex);
  }

  // If there are only two atoms to constrain, it does not make sense to constrain them.
  if (constrainedAtoms.size() <= 2)
    constrainedAtoms.clear();

  // Transfer the result for this subsystem to the ParametrizationData object
  data_.constrainedAtoms[subsystemIndex] = constrainedAtoms;
}

std::vector<int> ConstrainedAtomsIdentifier::getHeavyAtomTermini(const std::vector<unsigned>& indexMapForCurrentGraph,
                                                                 const Molassembler::Graph& graph,
                                                                 const Utils::AtomCollection& subsystem,
                                                                 const std::vector<int>& atomIndexMappingForCurrentSubsystem) {
  std::vector<int> heavyAtomTermini;
  // Loop over all atoms in the molecule
  for (auto atomIndex : graph.atoms()) {
    // Check if it is a heavy atom
    if (Utils::ElementInfo::Z(subsystem.getElement(indexMapForCurrentGraph[atomIndex])) > 1) {
      int numberOfAdjacentSaturatingHydrogens = 0;
      // Loop over all adjacent atoms
      for (auto adjacent : graph.adjacents(atomIndex)) {
        // Check if the adjacent atom is a saturating hydrogen
        if (atomIndexMappingForCurrentSubsystem[indexMapForCurrentGraph[adjacent]] == SwooseUtilities::indexForSaturatingAtoms) {
          numberOfAdjacentSaturatingHydrogens++;
        }
      }
      // Add heavy atom to the termini atoms if it has at least one saturating hydrogen bonded to it.
      if (numberOfAdjacentSaturatingHydrogens > 0)
        heavyAtomTermini.push_back(atomIndex);
    }
  }
  return heavyAtomTermini;
}

void ConstrainedAtomsIdentifier::updateConstrainedAtomsForOneMolecule(std::vector<int>& constrainedAtoms,
                                                                      const std::vector<int>& heavyAtomsToConstrain,
                                                                      const std::vector<unsigned>& indexMapForCurrentGraph,
                                                                      const Molassembler::Graph& graph, int subsystemIndex) {
  for (auto atomIndex : heavyAtomsToConstrain) {
    constrainedAtoms.push_back(indexMapForCurrentGraph[atomIndex]);
    // Loop over all adjacent atoms
    for (auto adjacent : graph.adjacents(atomIndex)) {
      // Check whether the adjacent atom is a saturating hydrogen
      if (data_.atomIndexMapping[subsystemIndex][indexMapForCurrentGraph[adjacent]] == SwooseUtilities::indexForSaturatingAtoms) {
        // Add it to the constrained atoms and then break, since we only want to constrain one bond per heavy atom
        constrainedAtoms.push_back(indexMapForCurrentGraph[adjacent]);
        break;
      }
    }
  }
}

} // namespace MMParametrization
} // namespace Scine
