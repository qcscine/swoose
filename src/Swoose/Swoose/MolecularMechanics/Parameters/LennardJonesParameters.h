/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_LENNARDJONESPARAMETERS_H
#define MOLECULARMECHANICS_LENNARDJONESPARAMETERS_H

#include "../Interactions/LennardJones.h"

namespace Scine {
namespace MolecularMechanics {
/**
 * @class LennardJonesParameters LennardJonesParameters.h
 * @brief Class containing the parameters for a van der Waals LJ-type interaction.
 */
class LennardJonesParameters {
 public:
  /**
   * @brief Constructor.
   * @param vdwRadius Unit: Angstrom
   * @param wellDepth Unit: kcal/mol
   */
  LennardJonesParameters(double vdwRadius, double wellDepth);

  /**
   * @brief Method returning the LennardJones analogon with the right unit conversions to be used in the calculation.
   */
  LennardJones toMMLennardJones(const LennardJonesParameters& otherLjParameters, double scalingFactor) const;

 private:
  double vdwRadius_; // Unit: Angstrom
  double wellDepth_; // Unit: kcal/mol
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_BONDPARAMETERS_H