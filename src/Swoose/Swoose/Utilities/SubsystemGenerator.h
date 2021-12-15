/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_SUBSYSTEMGENERATOR_H
#define SWOOSEUTILITIES_SUBSYSTEMGENERATOR_H

#include <list>
#include <memory>
#include <random>
#include <vector>

namespace Scine {

namespace Core {
class Log;
} // namespace Core

namespace Utils {
class Atom;
class AtomCollection;
class BondOrderCollection;
} // namespace Utils

namespace SwooseUtilities {
class FragmentAnalyzer;

class SubsystemGenerator {
 public:
  /**
   * @brief Constructor.
   *
   * @param fullStructure The molecular system's full structure.
   * @param bondOrders The bond order matrix of the full system.
   * @param fragmentAnalyzer A reference to a FragmentAnalyzer class.
   * @param bondOrderThreshold Threshold of which bond order to still consider as a bond.
   * @param maximumSubsystemSize Maximum size of a subsystem. If a large subsystem is generated, a warning is printed.
   * @param log The logger.
   * @param randomSeed The random seed that is set in the constructor. Default: 42
   * @param probabilityToDivideBond Probability to cleave a bond that is in principle cleavable. Default: 1.0
   */
  SubsystemGenerator(const Utils::AtomCollection& fullStructure, const Utils::BondOrderCollection& bondOrders,
                     FragmentAnalyzer& fragmentAnalyzer, double bondOrderThreshold, int maximumSubsystemSize,
                     Core::Log& log, int randomSeed = 42, double probabilityToDivideBond = 1.0);

  /**
   * @brief Generates a subsystem (fragment) from a given full system.
   *
   * @param centralAtomIndex The index of the central atom around which the subsystem will be constructed.
   * @param atomIndexMapping Vector of indices that correspond to the indices of the atoms
   *                         inside the given fragment in the full system. It is given as a reference and updated.
   *                         Therefore, it is usually provided to the function as an empty vector.
   * @param additionToRadius Variable extra amount that is added to the base subsystem radius and
   *                         is updated whenever necessary.
   */
  Utils::AtomCollection generateSubsystem(int centralAtomIndex, std::vector<int>& atomIndexMapping, double subsystemRadius);

 private:
  /*
   * @brief Tries to update the given subsystem generating a guess for a sensible (in terms of charge and multiplicity)
   *        subsystem structure. All parameters except for the atom index are passed by reference or const reference.
   *
   * @param subsystem The current subsystem.
   * @param centralAtom The atom in the center of the current subsystem (around which the fragment is build).
   * @param atomIndex The index of the atom at the center of the current subsystem (also subsystem index).
   * @param atomIndexMapping Vector of indices that correspond to the indices of the atoms
   *                         inside the given fragment in the full system. It is given as a reference and updated.
   *                         Therefore, it is usually provided to the function as an empty vector.
   * @param subsystemRadius The initial radius of the sphere defining the subsystem.
   * @param additionToRadius Variable extra amount that is added to the base subsystem radius and
   *                         is updated whenever necessary.
   * @param unsuccessful Boolean tracking whether the subsystem generation is still unsuccessful.
   */
  void tryGeneratingSensibleSubsystem(Utils::AtomCollection& subsystem, const Utils::Atom& centralAtom, int atomIndex,
                                      std::vector<int>& atomIndexMapping, const double& subsystemRadius,
                                      double& additionToRadius, bool& unsuccessful);
  /*
   * Vector of size N (number of atoms in the full system) with the following value for
   * each atom: the size of the molecular subgraph (individual molecule) it is in.
   */
  std::vector<int> subgraphSizes_;
  // The molecular structure of the whole system.
  const Utils::AtomCollection& fullStructure_;
  // Bond orders of the full system.
  const Utils::BondOrderCollection& bondOrders_;
  // The connectivity of the system. It is a vector of a list of neighbor atom indices for each atom.
  std::vector<std::list<int>> listsOfNeighbors_;
  // Fragment analyzer class reference
  FragmentAnalyzer& fragmentAnalyzer_;
  // Threshold of which bond order to still consider as a bond.
  double bondOrderThreshold_;
  // Maximum size of a subsystem. If a large subsystem is generated, a warning is printed.
  int maximumSubsystemSize_;
  // Random engine for non-deterministic fragment generation.
  std::shared_ptr<std::mt19937> randomEngine_;
  // Logger
  Core::Log& log_;
  // Probability to cleave a bond that is in principle cleavable.
  double probabilityToDivideBond_;
  // Maximum number of fragmentation attempts for one fragment before an exception is thrown
  static constexpr int maxNumOfAttemptsForOneSubsystem_ = 100;
};

} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_SUBSYSTEMGENERATOR_H
