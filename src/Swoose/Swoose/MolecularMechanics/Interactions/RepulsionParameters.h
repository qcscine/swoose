/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_REPULSIONPARAMETERS_H
#define MOLECULARMECHANICS_REPULSIONPARAMETERS_H

#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <memory>

namespace Scine {
namespace MolecularMechanics {

/**
 * @class RepulsionParameters RepulsionParameters.h
 * @brief Class handling the parameters needed for the repulsive non-bonded interaction term.
 */
class RepulsionParameters {
 public:
  /**
   * @brief Constructor.
   * @param R0 The matrix of D3 cutoff radii.
   * @param betaRepulsion The global beta parameter for the Pauli repulsion.
   */
  explicit RepulsionParameters(const Eigen::MatrixXd& R0, const double& betaRepulsion);
  /**
   * @brief Getter for the effective charge of a certain element.
   */
  double getEffectiveCharge(Utils::ElementType element) const;
  /**
   * @brief Getter for the D3 cutoff radius R0 for two atoms.
   */
  double getR0(int atom1Index, int atom2Index) const;
  /**
   * @brief Getter for beta parameter.
   */
  double getBetaRepulsion();

 private:
  int getValenceElectronNumbers(Utils::ElementType element) const;
  double getValenceElectronScalingFactor(Utils::ElementType element) const;

  // The matrix of D3 cutoff radii.
  const Eigen::MatrixXd& R0_;

  // The global beta parameter for the Pauli repulsion.
  double betaRepulsion_;

  // Taken from J. Chem. Theory Comput. 2014, 10, 4497-4514
  double valenceElectronScalingFactors_[94] = {
      2.35, 2.35, 1.7,  5.5,  0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 2.5,  3.0,  0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 3.0,
      3.0,  0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 3.0,  3.0,
      0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 0.63, 3.0,  3.0,  0.6,
      0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,
      0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6,  3.0,  3.0,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6};

  // Taken from wikipedia: https://en.wikipedia.org/wiki/List_of_elements_by_atomic_properties
  int valenceElectrons_[94] = {1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 2, 2, 2, 1,
                               2, 2, 2, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, 1, 2,
                               3, 4, 5, 6, 7, 8, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                               2, 2, 2, 2, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 2, 2, 2, 2, 2, 2};
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_REPULSIONPARAMETERS_H
