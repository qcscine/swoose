/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_DIHEDRALPARAMETERS_H
#define MOLECULARMECHANICS_DIHEDRALPARAMETERS_H

#include "../Interactions/Dihedral.h"

namespace Scine {
namespace MolecularMechanics {
/**
 * @class DihedralParameters DihedralParameters.h
 * @brief Class containing the parameters for an MM dihedral.
 */
class DihedralParameters {
 public:
  /**
   * @brief Constructor.
   * @param halfBarrierHeight The barrier height divided by a factor of 2. Unit: kcal/mol
   * @param phaseShift The phase shift angle in the torsional function. Unit: degrees
   * @param periodicity The periodicity of the torsional barrier (number of local minima) */
  DihedralParameters(double halfBarrierHeight, double phaseShift, int periodicity);

  /**
   * @brief Method returning the MMDihedral analogon with the right unit conversions to be used in the calculation.
   */
  Dihedral toMMDihedral() const;

  /**
   * @brief Tells if the dihedral has no contribution.
   */
  bool isZero() const;

  /** @brief Setter for the half barrier height */
  void setHalfBarrierHeight(const double& hbh);
  /** @brief Setter for the phase shift */
  void setPhaseShift(const double& ps);
  /** @brief Setter for the periodicity */
  void setPeriodicity(const int& p);
  /** @brief Getter for the half barrier height */
  double getHalfBarrierHeight() const;
  /** @brief Getter for the phase shift */
  double getPhaseShift() const;
  /** @brief Getter for the periodicity */
  int getPeriodicity() const;

 private:
  double halfBarrierHeight_; // unit: kcal/mol
  double phaseShift_;        // unit: degrees
  int periodicity_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_DIHEDRALPARAMETERS_H
