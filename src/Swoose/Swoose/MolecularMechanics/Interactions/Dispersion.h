/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_DISPERSION_H
#define MOLECULARMECHANICS_DISPERSION_H

#include <Utils/Math/AutomaticDifferentiation/Second1D.h>
#include <Utils/Typenames.h>

namespace Scine {
namespace MolecularMechanics {
/**
 * @class Dispersion Dispersion.h
 * @brief  Class treating a dispersion (D3-BJ) interaction, based solely on the bond length. (i.e. in 1 dimension)
 */
class Dispersion {
 public:
  /** @brief Constructor. */
  Dispersion(double scalingFactor, double c6);
  /** @brief Getter for scaling factor. */
  double getScalingFactor() const;
  /** @brief Getter for C6 coefficient. */
  double getC6() const;

 private:
  double scalingFactor_;
  double c6_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_DISPERSION_H
