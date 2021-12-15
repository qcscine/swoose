/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DihedralParameters.h"
#include <Utils/Constants.h>
#include <cmath>

namespace Scine {
namespace MolecularMechanics {

DihedralParameters::DihedralParameters(double halfBarrierHeight, double phaseShift, int periodicity)
  : halfBarrierHeight_(halfBarrierHeight), phaseShift_(phaseShift), periodicity_(periodicity) {
}

Dihedral DihedralParameters::toMMDihedral() const {
  return {halfBarrierHeight_ * Utils::Constants::hartree_per_kCalPerMol, periodicity_,
          phaseShift_ * Utils::Constants::rad_per_degree};
}

bool DihedralParameters::isZero() const {
  return std::abs(halfBarrierHeight_) < 1e-12;
}

void DihedralParameters::setHalfBarrierHeight(const double& hbh) {
  halfBarrierHeight_ = hbh;
}

void DihedralParameters::setPhaseShift(const double& ps) {
  phaseShift_ = ps;
}

void DihedralParameters::setPeriodicity(const int& p) {
  periodicity_ = p;
}

double DihedralParameters::getHalfBarrierHeight() const {
  return halfBarrierHeight_;
}

int DihedralParameters::getPeriodicity() const {
  return periodicity_;
}

double DihedralParameters::getPhaseShift() const {
  return phaseShift_;
}

} // namespace MolecularMechanics
} // namespace Scine
