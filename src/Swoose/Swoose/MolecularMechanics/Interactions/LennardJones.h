/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_LENNARDJONES_H
#define MOLECULARMECHANICS_LENNARDJONES_H

#include <Utils/Math/AutomaticDifferentiation/Second1D.h>

namespace Scine {
namespace MolecularMechanics {
/**
 * @class LennardJones LennardJones.h
 * @brief Class treating a non-bonded van der Waals interaction with a LJ potential,
 *        based solely on the bond length. (i.e. in 1 dimension)
 */
class LennardJones {
 public:
  /**
   * @brief Constructor with the a and b coefficients of the LJ potential according to the formulation
   *        at http://ambermd.org/vdwequation.pdf (visited 2021-02-09).
   */
  LennardJones(double aij, double bij);
  /**
   * @brief Constructor without arguments sets both parameters to zero and records that no parameters
   *        are available for this bond, hence throwing an exception if getInteraction() is called. This
   *        will not throw any error if the corresponding LJ term is disabled (e.g., in a QM/MM calculation).
   */
  LennardJones();
  /**
   * @brief Returns whether parameters are available for this LJ.
   */
  bool hasParameters() const;
  /**
   * @brief Evaluates the energy and its derivatives for one LJ term.
   */
  Utils::AutomaticDifferentiation::Second1D getInteraction(double distance) const;

 private:
  double aij_;
  double bij_;
  bool parametersAreAvailable_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_LENNARDJONES_H
