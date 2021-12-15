/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_FRAGMENTATIONHELPER_H
#define SWOOSEUTILITIES_FRAGMENTATIONHELPER_H

#include <Utils/Geometry/ElementTypes.h>
#include <deque>
#include <list>
#include <memory>
#include <random>
#include <vector>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace Utils {
class Settings;
class Atom;
class AtomCollection;
class BondOrderCollection;
} // namespace Utils

namespace SwooseUtilities {

// The index given to saturating atoms in the atomIndexMapping object
static constexpr int indexForSaturatingAtoms = -1;
// Minimum size of a subsystem
static constexpr int minimumSubsystemSize_ = 20;

namespace FragmentationHelper {
/**
 * @brief This function recursively searches for the right places to cut a fragment out of the system starting
 *        from the bond atomInside - atomOutside.
 * @param atomsToAdd AtomCollection with atoms that are later added to the raw subsystem
 * @param alreadyAddedAtoms Vector of atom indices in the full system for the atoms in atomsToAdd
 * @param isSaturatingAtom Deque of bools to decide whether an atom in atomsToAdd
 *        is an atom added for valence saturation.
 * @param atomInside Index of the atom inside of the current definition of the subsystem.
 * @param atomOutside Index of the atom outside of the current definition of the subsystem.
 * @param fullStructure The full system's molecular structure.
 * @param listsOfNeighbors The full system's connectivity as lists of neighbors.
 * @param probabilityToDivide Probability of dividing at a bond that is divisible in principle.
 * @param randomEngine A pointer to a Mersenne Twister pseudo-random generator.
 */
void addAtomsUpToReasonableCut(Utils::AtomCollection& atomsToAdd, std::vector<int>& alreadyAddedAtoms,
                               std::deque<bool>& isSaturatingAtom, int atomInside, int atomOutside,
                               const Utils::AtomCollection& fullStructure, const std::vector<std::list<int>>& listsOfNeighbors,
                               double probabilityToDivide, std::shared_ptr<std::mt19937> randomEngine);
/**
 * @brief This function adds a set of atoms to an existing molecular structure.
 * @param atomsToAdd The atoms to add.
 * @param subsystem The existing molecular structure.
 * @param isSaturatingAtom Deque of bools listing whether a given atom to add will be added for valence saturation.
 * @param addSaturatingAtoms Decides whether this function should add atoms placed for valence saturation.
 */
void addMoreAtomsToSubsystem(const Utils::AtomCollection& atomsToAdd, Utils::AtomCollection& subsystem,
                             const std::deque<bool>& isSaturatingAtom, bool addSaturatingAtoms);
/**
 * @brief Checks whether a fragment is too small or too large.
 *        If it is too small, the fragmentation is repeated with a larger radius,
 *        if it is too large, a warning is given.
 */
void checkSizeOfSubsystem(int size, double& additionToRadius, bool& unsuccessfulFragmentation, int fragmentIndex,
                          int fullStructureSize, int maximumSubsystemSize, Core::Log& log);
/**
 * @brief Updates the information of how the indices map from the subsystem to the full structure for a given subsystem.
 */
void updateInformationForIndexMapping(const Utils::AtomCollection& subsystem,
                                      const Utils::AtomCollection& fullStructure, std::vector<int>& atomIndexMapping);
/**
 * @brief Computes and returns a vector of size N (number of atoms in the full system) with the following value for
 *        each atom: the size of the molecular subgraph (individual molecule) it is in.
 */
std::vector<int> calculateSubgraphSizes(const Utils::AtomCollection& fullStructure, const Utils::BondOrderCollection& bondOrders);

/**
 * @brief Checks whether there is already an atom close to the given atom in the given molecular structure.
 * @param atom The atom to be evaluated.
 * @param existingAtoms The molecular structure.
 */
bool atomIsCloseToExistingAtom(const Utils::Atom& atom, const Utils::AtomCollection& existingAtoms);

/**
 * @brief Evaluates whether a bond connecting the two given atoms is a good place to cut.
 * @param atomInside Index of the atom that is already considered in the subsystem.
 * @param atomOutside Index of the atom that is not yet considered part of the subsystem.
 * @param fullStructure The full system's molecular structure.
 * @param listsOfNeighbors The connectivity of the system as lists of neighbors.
 * @param probabilityToDivide Probability of cleaving at a bond that is divisible in principle.
 * @param randomEngine A pointer to a Mersenne Twister pseudo-random generator.
 */
bool isDivisibleAtBond(int atomInside, int atomOutside, const Utils::AtomCollection& fullStructure,
                       const std::vector<std::list<int>>& listsOfNeighbors, double probabilityToDivide,
                       std::shared_ptr<std::mt19937> randomEngine);

} // namespace FragmentationHelper
} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_FRAGMENTATIONHELPER_H
