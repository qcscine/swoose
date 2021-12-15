/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DispersionTerm.h"
#include "DispersionEvaluator.h"
#include <Utils/Math/FullSecondDerivativeCollection.h>

namespace Scine {
namespace MolecularMechanics {
DispersionTerm::DispersionTerm(AtomIndex firstAtom, AtomIndex secondAtom, const Dispersion& dispersion,
                               std::shared_ptr<double> cutoffRadius)
  : firstAtom_(firstAtom), secondAtom_(secondAtom), dispersion_(dispersion), cutoffRadius_(cutoffRadius) {
}

DispersionTerm::~DispersionTerm() = default;

double DispersionTerm::evaluateDispersionTerm(const std::vector<Utils::Dftd3::Dftd3Atom>& structureOfDftd3Atoms,
                                              Utils::FullSecondDerivativeCollection& derivatives,
                                              std::shared_ptr<Utils::Dftd3::Dftd3> d3, Eigen::MatrixXd& R0) const {
  if (this->disabled_)
    return 0.0;

  auto rVector = structureOfDftd3Atoms[secondAtom_].getPosition() - structureOfDftd3Atoms[firstAtom_].getPosition();

  if (rVector.norm() > *cutoffRadius_)
    return 0.0;

  const double& factor = dispersion_.getScalingFactor();

  auto c6 = dispersion_.getC6();
  auto c8 = d3->calculateC8Coefficient(structureOfDftd3Atoms[firstAtom_], structureOfDftd3Atoms[secondAtom_], c6);
  auto r0 = sqrt(c8 / c6);
  R0(firstAtom_, secondAtom_) = r0;

  Utils::AutomaticDifferentiation::Second1D r(rVector.norm(), 1, 0);

  auto r2 = r * r;
  auto r3 = r * r2;
  auto r6 = r3 * r3;
  auto r8 = r6 * r2;

  auto e6 = c6 / (r6 + pow(DispersionEvaluator::getA1() * r0 + DispersionEvaluator::getA2(), 6));
  auto e8 = c8 / (r8 + pow(DispersionEvaluator::getA1() * r0 + DispersionEvaluator::getA2(), 8));

  Utils::AutomaticDifferentiation::Second1D energy = -factor * (e6 + DispersionEvaluator::getS8() * e8);

  auto v = Utils::AutomaticDifferentiation::get3Dfrom1D<Utils::DerivativeOrder::Two>(energy, rVector);

  derivatives.addDerivative(firstAtom_, secondAtom_, v);

  return energy.value();
}

int DispersionTerm::getFirstAtom() const {
  return firstAtom_;
}

int DispersionTerm::getSecondAtom() const {
  return secondAtom_;
}

Dispersion DispersionTerm::getDispersion() const {
  return dispersion_;
}

} // namespace MolecularMechanics
} // namespace Scine