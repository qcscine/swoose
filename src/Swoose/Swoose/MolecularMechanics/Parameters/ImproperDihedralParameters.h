/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_IMPROPERDIHEDRALPARAMETERS_H
#define MOLECULARMECHANICS_IMPROPERDIHEDRALPARAMETERS_H

#include "../Interactions/ImproperDihedral.h"

namespace Scine {
namespace MolecularMechanics {
/**
 * @class ImproperDihedralParameters ImproperDihedralParameters.h
 * @class Class containing the parameters for an MM improper dihedral.
 */
class ImproperDihedralParameters {
 public:
  /**
   * @brief Constructor.
   * @param forceConstant Unit: kcal/mol
   * @param equilibriumAngle Unit: degrees */
  ImproperDihedralParameters(double forceConstant, double equilibriumAngle);

  /**
   * @brief Method returning the MMImproperDihedral analogon to be used in the calculation.
   */
  ImproperDihedral toMMImproperDihedral() const;

  /** @brief Setter for the force constant */
  void setForceConstant(const double& fc);
  /** @brief Setter for the equilibrium angle */
  void setEquilibriumAngle(const double& eqAng);
  /** @brief Getter for the force constant */
  double getForceConstant() const;
  /** @brief Getter for the equilibrium angle */
  double getEquilibriumAngle() const;

 private:
  double forceConstant_;    // unit: kcal/mol
  double equilibriumAngle_; // unit: degrees
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_IMPROPERDIHEDRALPARAMETERS_H
