/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_ATOMICCHARGESASSEMBLER_H
#define MMPARAMETRIZATION_ATOMICCHARGESASSEMBLER_H

#include <vector>

namespace Scine {
namespace Core {
struct Log;
}

namespace MMParametrization {
struct ParametrizationData;

namespace AtomicChargesAssembler {

/**
 * @brief Assemble atomic charges of the full nanoscale system from the fragment calculations.
 * @param data The ParametrizationData object given as a reference.
 * @param log The logger.
 */
void assembleAtomicCharges(ParametrizationData& data, Core::Log& log);

/**
 * @brief Renormalizes the atomic charges in the given data object,
 *        so that the total charge of the system is correct.
 * @param data The ParametrizationData object given as a reference.
 */
void renormalizeAtomicCharges(ParametrizationData& data);

/**
 * @brief Updates the atomic charge for a given atom, which it takes from user-provided parameters already stored
 *        in the parametrization data.
 * @param atomicCharge The atomic charge to change as a reference.
 * @param atomIndex The atom index of the atom to which that atomic charge belongs.
 * @param data The ParametrizationData object given as a const reference.
 * @return Whether it was successful to obtain that charge (it is not if there were no user-provided parameters).
 */
bool updateChargeWithExistingParameters(double& atomicCharge, int atomIndex, const ParametrizationData& data);

} // namespace AtomicChargesAssembler
} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_ATOMICCHARGESASSEMBLER_H
