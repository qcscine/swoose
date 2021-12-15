/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_DISPERSIONEVALUATOR_H
#define MOLECULARMECHANICS_DISPERSIONEVALUATOR_H

#include "DispersionTerm.h"

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace Utils {
class AtomCollection;
class FullSecondDerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {

/**
 * @brief DispersionEvaluator DispersionEvaluator.h
 * @brief This class evaluates the overall energy and derivatives of the dispersion interactions.
 */
class DispersionEvaluator {
 public:
  /// @brief Constructor.
  explicit DispersionEvaluator(const Utils::AtomCollection& structure);
  /**
   * @brief This function evaluates and returns the energy for all dispersion interactions and updates the derivatives.
   * @param derivatives The derivatives.
   * @param d3 A pointer to an instance of the Dftd3 class.
   * @param R0 The matrix of D3 cutoff radii.
   * @return The dispersion energy.
   */
  double evaluate(Utils::FullSecondDerivativeCollection& derivatives, std::shared_ptr<Utils::Dftd3::Dftd3> d3,
                  Eigen::MatrixXd& R0);
  /**
   * @brief Sets a vector of instances of the DispersionTerm class.
   */
  void setDispersionTerms(std::vector<DispersionTerm>&& dispersionTerms);

  /// @brief Getter for the D3 parameter a1.
  static double getA1();
  /// @brief Getter for the D3 parameter s8.
  static double getS8();
  /// @brief Getter for the D3 parameter a2.
  static double getA2();

  /**
   * @brief Setter for the D3 parameters a1, s8 and a2.
   */
  void setD3Parameters(std::vector<double> d3Parameters);

 private:
  // friend class declaration
  friend class Qmmm::InteractionTermEliminator;
  const Utils::AtomCollection& structure_;
  std::vector<DispersionTerm> dispersions_;
  Eigen::MatrixXd dC6_dCN_;
  Eigen::MatrixXd dCN_drij_;
  // The three D3 parameters
  static double a1_;
  static double s8_;
  static double a2_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_DISPERSIONEVALUATOR_H
