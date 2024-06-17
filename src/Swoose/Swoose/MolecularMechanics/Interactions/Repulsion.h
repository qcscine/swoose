/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_REPULSION_H
#define MOLECULARMECHANICS_REPULSION_H

#include <Utils/Math/AutomaticDifferentiation/Second1D.h>

namespace Scine {
namespace MolecularMechanics {

/**
 * @class Repulsion Repulsion.h
 * @brief Class treating a repulsive non-bonded interaction, based solely on the bond length. (i.e. in 1 dimension)
 */

class Repulsion {
 public:
  Repulsion(double scalingFactor);

  Utils::AutomaticDifferentiation::Second1D getInteraction(const double& bondLength, const double& effectiveChargeA,
                                                           const double& effectiveChargeB, const double& betaRepulsion,
                                                           const double& R0) const;

  double scalingFactor_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_REPULSION_H
