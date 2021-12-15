/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_PARAMETEROPTIMIZER_H
#define MMPARAMETRIZATION_PARAMETEROPTIMIZER_H

#include <memory>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace Utils {
class Settings;
} // namespace Utils

namespace MMParametrization {
struct ParametrizationData;

/**
 * @class ParameterOptimizer ParameterOptimizer.h
 * @brief The class optimizes the force constants of the MM model in a least squares manner.
 */
class ParameterOptimizer {
 public:
  /**
   * @brief Constructor.
   */
  ParameterOptimizer(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings, Core::Log& log);
  /**
   * @brief This function optimizes the force constants of the MM model, which part of the MM parameters
   *        stored in the "parameters" member of "data_".
   */
  void optimizeParameters();

 private:
  // Implementation of the optimization function
  void optimizeParametersImpl();
  // Implementation of the optimization of dihedral parameters
  void optimizeDihedralParametersImpl();
  // Implementation of the optimization of angle parameters
  void optimizeAngleParametersImpl();
  // Implementation of the optimization of bond parameters
  void optimizeBondParametersImpl();
  // Implementation of the optimization of improper dihedral parameters
  void optimizeImproperDihedralParametersImpl();
  // The data used within all MM parametrization classes
  ParametrizationData& data_;
  // The settings
  std::shared_ptr<Utils::Settings> settings_;
  // The logger.
  Core::Log& log_;
  // Maximum number of function evaluations employed for the least squares partial Hessian fit
  static constexpr int maximumNumberOfFunctionEvaluationsInOptimization_ = 7;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_PARAMETEROPTIMIZER_H
