/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "HydrogenBondTerm.h"
#include "../MMExceptions.h"
#include "AngleTerm.h"
#include "HydrogenBondParameters.h"
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/FullSecondDerivativeCollection.h>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace MolecularMechanics {

HydrogenBondTerm::HydrogenBondTerm(AtomIndex donorAtom, AtomIndex hydrogenAtom, AtomIndex acceptorAtom,
                                   const HydrogenBond& hydrogenBond)
  : donorAtom_(donorAtom), hydrogenAtom_(hydrogenAtom), acceptorAtom_(acceptorAtom), hydrogenBond_(hydrogenBond) {
}

HydrogenBondTerm::~HydrogenBondTerm() = default;

double HydrogenBondTerm::evaluateHydrogenBondTerm(const Utils::AtomCollection& structure,
                                                  Utils::AtomicSecondDerivativeCollection& derivatives,
                                                  const std::vector<double>& atomicCharges) const {
  if (this->disabled_)
    return 0.0;

  const auto& elements = structure.getElements();
  auto positions = structure.getPositions();
  Eigen::Vector3d a(positions.row(donorAtom_) - positions.row(hydrogenAtom_));
  Eigen::Vector3d b(positions.row(acceptorAtom_) - positions.row(hydrogenAtom_));

  double angle = acos((a.dot(b) / (a.norm() * b.norm())));

  auto R = positions.row(acceptorAtom_) - positions.row(donorAtom_);

  double chargeDonor = atomicCharges[donorAtom_];
  double chargeAcceptor = atomicCharges[acceptorAtom_];

  double constantDonor = HydrogenBondParameters::getInteractionStrengthConstants(elements[donorAtom_]);
  double constantAcceptor = HydrogenBondParameters::getInteractionStrengthConstants(elements[acceptorAtom_]);

  Utils::AutomaticDifferentiation::Second1D resultDistanceVariable = hydrogenBond_.getInteractionDistanceVariable(
      R.norm(), angle, chargeDonor, chargeAcceptor, constantDonor, constantAcceptor);
  Utils::AutomaticDifferentiation::Second1D resultAngleVariable =
      hydrogenBond_.getInteractionAngleVariable(R.norm(), angle, chargeDonor, chargeAcceptor, constantDonor, constantAcceptor);

  if (std::abs(resultDistanceVariable.value() - resultAngleVariable.value()) > 1e-12)
    throw TwoResultsForHydrogenBondsAreNotEqualException(
        "The two results for the hydrogen bond interaction should be equal, but are not!");

  // Add derivative contribution for the variable distance (R)
  auto v = Utils::AutomaticDifferentiation::get3Dfrom1D<Utils::DerivativeOrder::Two>(resultDistanceVariable, R);
  derivatives[acceptorAtom_] += v;
  derivatives[donorAtom_] += v.opposite();

  // Add derivative contribution for the variable angle
  Utils::AutomaticDifferentiation::Second3D angleContributionAtomDonor;
  Utils::AutomaticDifferentiation::Second3D angleContributionAtomHydrogen;
  Utils::AutomaticDifferentiation::Second3D angleContributionAtomAcceptor;
  if (angle < Utils::Constants::pi - AngleTerm::singularityCriterion_ && angle > AngleTerm::singularityCriterion_) {
    AngleTerm::calculateDerivativesWithNormalFormula(resultAngleVariable, a, b, angleContributionAtomDonor,
                                                     angleContributionAtomHydrogen, angleContributionAtomAcceptor);
  }
  else {
    AngleTerm::calculateDerivativesForCriticalAngles(resultAngleVariable, a, b, angleContributionAtomDonor,
                                                     angleContributionAtomHydrogen, angleContributionAtomAcceptor);
  }

  derivatives[donorAtom_] += angleContributionAtomDonor;
  derivatives[hydrogenAtom_] += angleContributionAtomHydrogen;
  derivatives[acceptorAtom_] += angleContributionAtomAcceptor;

  return resultDistanceVariable.value();
}

int HydrogenBondTerm::getDonorAtom() const {
  return donorAtom_;
}

int HydrogenBondTerm::getHydrogenAtom() const {
  return hydrogenAtom_;
}

int HydrogenBondTerm::getAcceptorAtom() const {
  return acceptorAtom_;
}

} // namespace MolecularMechanics
} // namespace Scine