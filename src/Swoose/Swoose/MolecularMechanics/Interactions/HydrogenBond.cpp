/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "HydrogenBond.h"
#include "Utils/Constants.h"
#include <cmath>

namespace Scine {
namespace MolecularMechanics {

static constexpr double dampingCutoff = 4.0 * Utils::Constants::bohr_per_angstrom;
static constexpr double kq1 = 10.0;
static constexpr double kq2 = 5.0;
static constexpr double interactionStrengthScalingFactor = 1.0;

HydrogenBond::HydrogenBond() {
}

Utils::AutomaticDifferentiation::Second1D
HydrogenBond::getInteractionDistanceVariable(double distance, double angle, double chargeDonor, double chargeAcceptor,
                                             double constantDonor, double constantAcceptor) const {
  Utils::AutomaticDifferentiation::Second1D r(distance, 1.0, 0.0);
  auto r2 = r * r;
  auto r3 = r * r2;
  auto r6 = r3 * r3;
  auto r12 = r6 * r6;

  double dampingAngle = pow(0.5 * (cos(angle + Utils::Constants::pi) + 1), 6);

  Utils::AutomaticDifferentiation::Second1D dampingDistance = 1.0 / (1.0 + (r12 / pow(dampingCutoff, 12)));

  double interactionStrengthDonor =
      interactionStrengthScalingFactor * constantDonor * exp(-kq1 * chargeDonor) / (exp(-kq1 * chargeDonor) + kq2);
  double interactionStrengthAcceptor = interactionStrengthScalingFactor * constantAcceptor *
                                       exp(-kq1 * chargeAcceptor) / (exp(-kq1 * chargeAcceptor) + kq2);
  double interactionStrengthSum = interactionStrengthDonor + interactionStrengthAcceptor;

  return -1.0 * dampingAngle * dampingDistance * interactionStrengthSum * (1.0 / r3);
}

Utils::AutomaticDifferentiation::Second1D
HydrogenBond::getInteractionAngleVariable(double distance, double angle, double chargeDonor, double chargeAcceptor,
                                          double constantDonor, double constantAcceptor) const {
  Utils::AutomaticDifferentiation::Second1D theta(angle, 1.0, 0.0);

  Utils::AutomaticDifferentiation::Second1D dampingAngleFactor = 0.5 * (cos(theta + Utils::Constants::pi) + 1);
  auto dampingAngleFactor2 = dampingAngleFactor * dampingAngleFactor;
  auto dampingAngleFactor3 = dampingAngleFactor * dampingAngleFactor2;
  auto dampingAngle = dampingAngleFactor3 * dampingAngleFactor3;

  double r3 = distance * distance * distance;
  double r12 = pow(distance, 12);
  double dampingDistance = 1.0 / (1.0 + (r12 / pow(dampingCutoff, 12)));

  double interactionStrengthDonor =
      interactionStrengthScalingFactor * constantDonor * exp(-kq1 * chargeDonor) / (exp(-kq1 * chargeDonor) + kq2);
  double interactionStrengthAcceptor = interactionStrengthScalingFactor * constantAcceptor *
                                       exp(-kq1 * chargeAcceptor) / (exp(-kq1 * chargeAcceptor) + kq2);
  double interactionStrengthSum = interactionStrengthDonor + interactionStrengthAcceptor;

  return -1.0 * dampingAngle * dampingDistance * interactionStrengthSum * (1.0 / r3);
}

} // namespace MolecularMechanics
} // namespace Scine
