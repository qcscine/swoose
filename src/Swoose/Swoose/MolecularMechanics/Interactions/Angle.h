/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_ANGLE_H
#define MOLECULARMECHANICS_ANGLE_H

#include <Utils/Math/AutomaticDifferentiation/Second1D.h>

namespace Scine {
namespace MolecularMechanics {
/**
 * @class Angle Angle.h
 * @brief Class treating an angle interaction,based solely on the angle (in rad), i.e. in 1 dimension.
 */
class Angle {
 public:
  /** @brief Constructor from an equilibrium angle and a force constant. */
  Angle(double equilibriumAngle, double forceConstant);
  /**
   * @brief Constructor without arguments sets both parameters to zero and records that no parameters
   *        are available for this angle, hence throwing an exception if getInteraction() is called. This
   *        will not throw any error if the corresponding angle term is disabled (e.g., in a QM/MM calculation).
   */
  Angle();
  /**
   * @brief Returns whether parameters are available for this angle.
   */
  bool hasParameters() const;
  /** @brief Calculates the energy with derivatives of one angle term based on the angle. */
  Utils::AutomaticDifferentiation::Second1D getInteraction(double angle) const;

 private:
  double equilibriumAngle_;
  double forceConstant_;
  bool parametersAreAvailable_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_ANGLE_H
