/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_CONSTRAINEDATOMSIDENTIFIER_H
#define MMPARAMETRIZATION_CONSTRAINEDATOMSIDENTIFIER_H

#include <vector>

namespace Scine {

namespace Molassembler {
class Graph;
} // namespace Molassembler

namespace Utils {
class AtomCollection;
}

namespace MMParametrization {
struct ParametrizationData;

class ConstrainedAtomsIdentifier {
 public:
  /**
   * @brief Constructor.
   */
  explicit ConstrainedAtomsIdentifier(ParametrizationData& data);
  /**
   * @brief Updates the information about the constrained atoms in the geometry optimizations in the
   *        ParametrizationData object.
   * @param subsystem The given molecular subsystem.
   * @param subsystemIndex The index of the given molecular subsystem.
   */
  void updateInformationAboutConstrainedAtoms(const Utils::AtomCollection& subsystem, int subsystemIndex);

 private:
  /*
   * @brief Returns a vector of atom indices corresponding to all non-hydrogen atoms that are next to
   *        a bond, which was cleft during the fragmentation.
   */
  std::vector<int> getHeavyAtomTermini(const std::vector<unsigned>& indexMapForCurrentGraph,
                                       const Molassembler::Graph& graph, const Utils::AtomCollection& subsystem,
                                       const std::vector<int>& atomIndexMappingForCurrentSubsystem);
  /*
   * @brief Updates the current vector of atom indices of the constrained atoms, which is given as the parameter
   *        "constrainedAtoms" (passed by reference). This update corresponds to one specific molecule within the
   *        a given fragment with the index subsystemIndex.
   */
  void updateConstrainedAtomsForOneMolecule(std::vector<int>& constrainedAtoms, const std::vector<int>& heavyAtomsToConstrain,
                                            const std::vector<unsigned>& indexMapForCurrentGraph,
                                            const Molassembler::Graph& graph, int subsystemIndex);
  // The data used within all MM parametrization classes
  ParametrizationData& data_;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_CONSTRAINEDATOMSIDENTIFIER_H
