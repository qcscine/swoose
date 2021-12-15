/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "LennardJonesParameters.h"
#include <Utils/Constants.h>
#include <cmath>

namespace Scine {
namespace MolecularMechanics {

LennardJonesParameters::LennardJonesParameters(double vdwRadius, double wellDepth)
  : vdwRadius_(vdwRadius), wellDepth_(wellDepth) {
}

/*
 * Calculation according to: http://ambermd.org/vdwequation.pdf (visited 2021-02-09)
 */
LennardJones LennardJonesParameters::toMMLennardJones(const LennardJonesParameters& otherLjParameters,
                                                      double scalingFactor) const {
  double rij = (vdwRadius_ + otherLjParameters.vdwRadius_) * Utils::Constants::bohr_per_angstrom;
  double eij = std::sqrt(wellDepth_ * otherLjParameters.wellDepth_) * Utils::Constants::hartree_per_kCalPerMol * scalingFactor;

  double rij3 = rij * rij * rij;
  double rij6 = rij3 * rij3;
  double rij12 = rij6 * rij6;

  double aCoeff = eij * rij12;
  double bCoeff = 2 * eij * rij6;

  return {aCoeff, bCoeff};
}

} // namespace MolecularMechanics
} // namespace Scine