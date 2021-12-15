/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <vector>

#ifndef MMPARAMETRIZATION_FRAGMENTDATADISTRIBUTOR_H
#  define MMPARAMETRIZATION_FRAGMENTDATADISTRIBUTOR_H

namespace Scine {
namespace MMParametrization {
struct ParametrizationData;

/**
 * @class FragmentDataDistributor FragmentDataDistributor.h
 * @brief This class handles the distribution of the available data for the parameter optimization.
 *        This mainly includes the task of providing a list of candidate fragments from which to take reference data
 *        for the parameters associated with a given atom.
 */
class FragmentDataDistributor {
 public:
  /**
   * @brief Constructor.
   * @param data The parametrization data object used throughout the whole parametrization process.
   */
  explicit FragmentDataDistributor(ParametrizationData& data);
  /**
   * @brief Checks whether the reference data that has already been collected is sufficient to
   *        perform the parametrization. Returns true, if this is the case.
   * @param refineConnectivity Decides whether reference bond orders are calculated during the parametrization process.
   * @param failedCalculations A vector containing integers that denote which of the calculations failed for a given
   *                           fragment: Multiple of 2 -> Hessian failed, Multiple of 3 -> Atomic charges failed,
   *                           Multiple of 5 -> Bond orders failed. Default: empty list -> this check will not be
   *                           performed.
   * @throws std::runtime_error Throws if the parametrization cannot be completed anymore
   *                            because too many calculations already failed.
   */
  bool referenceDataIsSufficient(bool refineConnectivity, std::vector<int> failedCalculations = {}) const;
  /**
   * @brief Returns the indices of candidate fragments to get the data from
   *        for the given fragment index 'fragmentIndex'.
   */
  std::vector<int> getCandidateFragments(int fragmentIndex) const;
  /**
   * @brief Updates a given list of candidate fragments to get the data from for the given
   *        fragment index 'fragmentIndex'.
   */
  void updateCandidateFragments(int fragmentIndex, std::vector<int>& listOfCandidates) const;
  /**
   * @brief Updates a given list of candidate fragments with third shell neighboring atoms
   *        (neighbors of neighbors of neighbors) to get the data from for the given
   *        fragment index 'fragmentIndex'. This can be necessary for terminal atoms, because their parameters
   *        may otherwise be depending only on a small number of candidate fragments, of which one must lead to
   *        convergenced reference calculations.
   */
  void updateCandidateFragmentsWithThirdShellNeighbors(int fragmentIndex, std::vector<int>& listOfCandidates) const;

 private:
  /*
   * @brief Returns whether enough data has been calculated already based on the given candidate fragments.
   * @throws std::runtime_error Throws if the parametrization cannot be completed anymore
   *                            because too many calculations already failed.
   */
  bool fragmentIsCoveredByData(const std::vector<int>& candidates, bool refineConnectivity,
                               const std::vector<int>& failedCalculations) const;
  /*
   * Returns whether the fragment with the given candidates can not be covered anymore because too many calculations
   * have failed already.
   */
  bool fragmentIsHopeless(const std::vector<int>& candidates, const std::vector<int>& failedCalculations) const;
  // The parametrization data
  ParametrizationData& data_;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_FRAGMENTDATADISTRIBUTOR_H
