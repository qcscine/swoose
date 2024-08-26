/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_REPARAMETRIZATIONHELPER_H
#define MMPARAMETRIZATION_REPARAMETRIZATIONHELPER_H

#include <string>
#include <unordered_set>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace MMParametrization {
struct ParametrizationData;

/**
 * @class ReparametrizationHelper ReparametrizationHelper.h
 * @brief This class provides functionalities for the re-parametrization of a system, which
 *        has only slightly been modified. The most of the parameters can be re-used without
 *        new reference data generation or parameter optimization.
 */
class ReparametrizationHelper {
 public:
  /**
   * @brief Constructor.
   */
  ReparametrizationHelper(ParametrizationData& data, Core::Log& log);

  /**
   * @brief Parses already provided parameters and stores them into the ParametrizationData object.
   * @param parameterFile The path to the parameter file of the provided parameters.
   */
  void parseProvidedParameters(const std::string& parameterFile);

  /**
   * @brief Removes all of the topological elements from the topology in the ParametrizationData object,
   *        of which the parameters are already covered by the provided ones.
   */
  void manipulateTopology();

  /**
   * @brief Returns whether reference data has to be calculated for a fragment with the given fragment index
   *        in the case of a model re-parametrization.
   * @param fragmentIndex Index of a fragment.
   */
  bool isRelevantFragment(int fragmentIndex) const;

 private:
  // The data used within all MM parametrization classes
  ParametrizationData& data_;
  // The logger.
  Core::Log& log_;
  /*
   * Container for all atoms whose fragments will be relevant during the re-parametrization process.
   * It is filled when executing the manipulateTopology() function and used by the SuperfluousFragmentIdentifier
   * class to identify for which fragments it is not necessary to calculate reference data.
   */
  std::unordered_set<int> relevantFragments_;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_REPARAMETRIZATIONHELPER_H
