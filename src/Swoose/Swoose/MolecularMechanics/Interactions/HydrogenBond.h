/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_HYDROGENBOND_H
#define MOLECULARMECHANICS_HYDROGENBOND_H

#include <Utils/Math/AutomaticDifferentiation/Second1D.h>

namespace Scine {
namespace MolecularMechanics {
/**
 * @class HydrogenBond HydrogenBond.h
 * @brief Class calculating the energy and derivatives
 *        for a hydrogen bond based solely on the distance or on the angle, i.e. in 1 dimension, respectively.
 */
class HydrogenBond {
 public:
  /** @brief Default Constructor. */
  HydrogenBond();

  /**
   * @brief Calculate energy contribution from the distance with derivatives.
   */
  Utils::AutomaticDifferentiation::Second1D getInteractionDistanceVariable(double distance, double angle, double chargeDonor,
                                                                           double chargeAcceptor, double constantDonor,
                                                                           double constantAcceptor) const;
  /**
   * @brief Calculate energy contribution from
   * the angle with derivatives.
   */
  Utils::AutomaticDifferentiation::Second1D getInteractionAngleVariable(double distance, double angle, double chargeDonor,
                                                                        double chargeAcceptor, double constantDonor,
                                                                        double constantAcceptor) const;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_HYDROGENBOND_H
