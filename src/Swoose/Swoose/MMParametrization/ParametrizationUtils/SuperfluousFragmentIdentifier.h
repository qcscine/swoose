/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_SUPERFLUOUSFRAGMENTIDENTIFIER_H
#define MMPARAMETRIZATION_SUPERFLUOUSFRAGMENTIDENTIFIER_H

#include <memory>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace MMParametrization {
struct ParametrizationData;
class ReparametrizationHelper;

namespace SuperfluousFragmentIdentifier {

/**
 * @brief Identifies the fragments for which reference data does not have to be calculated. The reason is that
 *        it is identical to another fragment belonging to a neighboring atom, and therefore, the other fragment
 *        will be utilized as a source of the reference data anyways.
 *
 * @param data The parametrization data object.
 * @param log The logger.
 * @param reparametrizationHelper A pointer to the ReparametrizationHelper class, which provides the information
 *                                whether a given fragment is necessary to compute based on already provided
 *                                parameters. If no such information is given, this pointer is a nullptr.
 */
void identifySuperfluousFragments(ParametrizationData& data, Core::Log& log,
                                  std::shared_ptr<ReparametrizationHelper> reparametrizationHelper);

} // namespace SuperfluousFragmentIdentifier
} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_SUPERFLUOUSFRAGMENTIDENTIFIER_H
