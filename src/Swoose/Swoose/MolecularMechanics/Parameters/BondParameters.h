/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_BONDPARAMETERS_H
#define MOLECULARMECHANICS_BONDPARAMETERS_H

#include "../Interactions/Bond.h"

namespace Scine {
namespace MolecularMechanics {
/**
 * @class BondParameters BondParameters.h
 * @brief Class containing the parameters for an MM bond.
 */
class BondParameters {
 public:
  /**
   * @brief Constructor.
   * @param forceConstant Unit: kcal/(mol*(Angstrom^2))
   * @param equilibriumBondLength Unit: Angstrom */
  BondParameters(double forceConstant, double equilibriumBondLength);

  /**
   * @brief Method returning the MMBond analogon with the right unit conversions to be used in the calculation.
   */
  Bond toMMBond() const;

  /** @brief Setter for the force constant */
  void setForceConstant(const double& fc);
  /** @brief Setter for the equilibrium bond length */
  void setEquilibriumBondLength(const double& eqBoLen);
  /** @brief Getter for the force constant */
  double getForceConstant() const;
  /** @brief Getter for the force constantequilibrium bond length */
  double getEquilibriumBondLength() const;

 private:
  double forceConstant_;         // Unit: kcal/(mol*(Angstrom^2))
  double equilibriumBondLength_; // Unit: Angstrom
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_BONDPARAMETERS_H
