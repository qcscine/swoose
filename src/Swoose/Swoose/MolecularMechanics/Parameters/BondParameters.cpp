/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BondParameters.h"
#include <Utils/Constants.h>

namespace Scine {
namespace MolecularMechanics {
namespace Constants {
static constexpr double forceConstantConversionFactor =
    Utils::Constants::hartree_per_kCalPerMol / (Utils::Constants::bohr_per_angstrom * Utils::Constants::bohr_per_angstrom);
} // namespace Constants

BondParameters::BondParameters(double forceConstant, double equilibriumBondLength)
  : forceConstant_(forceConstant), equilibriumBondLength_(equilibriumBondLength) {
}

Bond BondParameters::toMMBond() const {
  return {equilibriumBondLength_ * Utils::Constants::bohr_per_angstrom,
          forceConstant_ * Constants::forceConstantConversionFactor};
}

void BondParameters::setForceConstant(const double& fc) {
  forceConstant_ = fc;
}

void BondParameters::setEquilibriumBondLength(const double& eqBoLen) {
  equilibriumBondLength_ = eqBoLen;
}

double BondParameters::getForceConstant() const {
  return forceConstant_;
}

double BondParameters::getEquilibriumBondLength() const {
  return equilibriumBondLength_;
}

} // namespace MolecularMechanics
} // namespace Scine
