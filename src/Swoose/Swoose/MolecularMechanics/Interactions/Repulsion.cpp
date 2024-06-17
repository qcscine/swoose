/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Repulsion.h"

namespace Scine {
namespace MolecularMechanics {

Repulsion::Repulsion(double scalingFactor) : scalingFactor_(scalingFactor) {
}

Utils::AutomaticDifferentiation::Second1D Repulsion::getInteraction(const double& bondLength, const double& effectiveChargeA,
                                                                    const double& effectiveChargeB,
                                                                    const double& betaRepulsion, const double& R0) const {
  Utils::AutomaticDifferentiation::Second1D dist(bondLength, 1, 0);
  return scalingFactor_ * effectiveChargeA * effectiveChargeB * (1 / dist) * exp(-betaRepulsion * dist / R0);
}

} // namespace MolecularMechanics
} // namespace Scine
