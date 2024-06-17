/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Bond.h"

namespace Scine {
namespace MolecularMechanics {

Bond::Bond(double equilibriumDistance, double forceConstant)
  : equilibriumDistance_(equilibriumDistance), forceConstant_(forceConstant), parametersAreAvailable_(true) {
}

Bond::Bond() : equilibriumDistance_(0.0), forceConstant_(0.0), parametersAreAvailable_(false) {
}

bool Bond::hasParameters() const {
  return parametersAreAvailable_;
}

Utils::AutomaticDifferentiation::Second1D Bond::getInteraction(double bondLength) const {
  Utils::AutomaticDifferentiation::Second1D dist(bondLength - equilibriumDistance_, 1, 0);
  return 0.5 * forceConstant_ * dist * dist;
}

} // namespace MolecularMechanics
} // namespace Scine
