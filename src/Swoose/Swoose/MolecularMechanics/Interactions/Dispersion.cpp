/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Dispersion.h"

namespace Scine {
namespace MolecularMechanics {

Dispersion::Dispersion(double scalingFactor, double c6) : scalingFactor_(scalingFactor), c6_(c6) {
}

double Dispersion::getScalingFactor() const {
  return scalingFactor_;
}

double Dispersion::getC6() const {
  return c6_;
}

} // namespace MolecularMechanics
} // namespace Scine
