/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_ANGLEPARAMETERS_H
#define MOLECULARMECHANICS_ANGLEPARAMETERS_H

#include "../Interactions/Angle.h"

namespace Scine {
namespace MolecularMechanics {
/**
 * @class AngleParameters AngleParameters.h
 * @brief Class containing the parameters for an MM angle.
 */
class AngleParameters {
 public:
  /**
   * @brief Constructor.
   * @param forceConstant Unit: kcal/(mol*(rad^2))
   * @param equilibriumAngle Unit: degrees */
  AngleParameters(double forceConstant, double equilibriumAngle);

  /**
   * @brief Method returning the MMAngle analogon with the right unit conversion to be used in the calculation.
   */
  Angle toMMAngle() const;

  /** @brief Setter for the force constant. */
  void setForceConstant(const double& fc);
  /** @brief Setter for the equilibrium angle. */
  void setEquilibriumAngle(const double& eqAng);
  /** @brief Getter for the force constant. */
  double getForceConstant() const;
  /** @brief Getter for the equilibrium angle. */
  double getEquilibriumAngle() const;

 private:
  double forceConstant_;    // Unit: kcal/mol/(rad^2)
  double equilibriumAngle_; // Unit: degrees
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_ANGLEPARAMETERS_H
