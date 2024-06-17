/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_ELECTROSTATIC_H
#define MOLECULARMECHANICS_ELECTROSTATIC_H

#include <Utils/Math/AutomaticDifferentiation/Second1D.h>

namespace Scine {
namespace MolecularMechanics {
/**
 * @class Electrostatic Electrostatic.h
 * @brief Class treating non-bonded electrostatic interaction, based solely on the bond length. (i.e. in 1 dimension)
 */
class Electrostatic {
 public:
  /**
   * @brief Constructor.
   * @param scalingFactor Scaling factor for the electrostatic interaction.
   */
  explicit Electrostatic(double scalingFactor);

  /**
   * @brief This function calculates the electrostatic energy and its derivatives
   *        for two point charges at a certain distance.
   */
  Utils::AutomaticDifferentiation::Second1D getInteraction(double distance, double charge1, double charge2) const;

 private:
  double scalingFactor_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_ELECTROSTATIC_H
