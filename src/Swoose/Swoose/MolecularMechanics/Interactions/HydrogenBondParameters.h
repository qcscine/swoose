/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_HYDROGENBONDPARAMETERS_H
#define MOLECULARMECHANICS_HYDROGENBONDPARAMETERS_H

#include "../MMExceptions.h"
#include <Utils/Geometry/ElementTypes.h>

namespace Scine {
namespace MolecularMechanics {

class HydrogenBondParameters {
 public:
  static double getInteractionStrengthConstants(const Utils::ElementType& elementType);
};

// Fitted by manual inspection to several small hydrogen bonded dimer interaction energies
inline double HydrogenBondParameters::getInteractionStrengthConstants(const Utils::ElementType& elementType) {
  double constant;
  switch (elementType) {
    case Utils::ElementType::N:
      constant = 0.6;
      break;
    case Utils::ElementType::O:
      constant = 0.7;
      break;
    case Utils::ElementType::F:
      constant = 3.2;
      break;
    case Utils::ElementType::Cl:
      constant = 4.2;
      break;
    default:
      throw ElementTypeIsNotAllowedForHydrogenBondsException(
          "The given element type can not be used for hydrogen bonds.");
  }
  return constant;
}

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_HYDROGENBONDPARAMETERS_H
