/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LennardJonesParameters.h"
#include <Utils/Constants.h>
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <cmath>

namespace Scine {
namespace MolecularMechanics {

LennardJonesParameters::LennardJonesParameters(double vdwRadius, double wellDepth)
  : vdwRadius_(vdwRadius * Utils::Constants::bohr_per_angstrom),
    wellDepth_(wellDepth * Utils::Constants::hartree_per_kCalPerMol) {
}

double LennardJonesParameters::getVdwRadius() const {
  return vdwRadius_;
}

double LennardJonesParameters::getWellDepth() const {
  return wellDepth_;
}
unsigned int LennardJonesParameters::getAtomTypeIndex() const {
  return atomTypeIndex_;
}
void LennardJonesParameters::setAtomTypeIndex(unsigned int index) {
  atomTypeIndex_ = index;
}

LennardJonesPairParameters::LennardJonesPairParameters(const LennardJonesParameters& first,
                                                       const LennardJonesParameters& second) {
  /*
   * Calculation according to: http://ambermd.org/vdwequation.pdf (visited 2021-02-09)
   */
  double rij = first.getVdwRadius() + second.getVdwRadius();
  double eij = std::sqrt(first.getWellDepth() * second.getWellDepth());

  const double rij3 = rij * rij * rij;
  const double rij6 = rij3 * rij3;
  const double rij12 = rij6 * rij6;

  aCoeff_ = eij * rij12;
  bCoeff_ = 2 * eij * rij6;
}
const double& LennardJonesPairParameters::getACoeff() {
  return aCoeff_;
}
const double& LennardJonesPairParameters::getBCoeff() {
  return bCoeff_;
}
Utils::AutomaticDifferentiation::Second1D LennardJonesPairParameters::getInteraction(const double& distance,
                                                                                     const double& scaling) const {
  auto dist2 = distance * distance;
  auto dist3 = dist2 * distance;
  auto dist6 = dist3 * dist3;

  auto dist6inv = Utils::AutomaticDifferentiation::getFromFull<Utils::DerivativeOrder::Two>(
      1.0 / dist6, -6.0 / (dist6 * distance), 42.0 / (dist6 * dist2));
  auto dist12inv = dist6inv * dist6inv;

  return (aCoeff_ * dist12inv - bCoeff_ * dist6inv) * scaling;
}
} // namespace MolecularMechanics
} // namespace Scine
