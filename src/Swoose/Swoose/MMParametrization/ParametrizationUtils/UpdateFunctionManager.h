/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_UPDATEFUNCTIONMANAGER_H
#define MMPARAMETRIZATION_UPDATEFUNCTIONMANAGER_H

#include <Swoose/MolecularMechanics/SFAM/SfamMolecularMechanicsCalculator.h>
#include <Utils/Optimizer/LeastSquares/UpdateFunctionManagerBase.h>
#include <memory>

namespace Scine {

namespace Utils {
class Settings;
} // namespace Utils

namespace MMParametrization {
struct ParametrizationData;

/**
 * @brief The types of parameters are assigned a bit each to switch them on
 *        if they are present in the parameters being optimized.
 */
enum class ParameterToOptimize { Bond = (1 << 0), Angle = (1 << 1), Dihedral = (1 << 2), ImproperDihedral = (1 << 3) };
/**
 * @brief Returns a ParameterToOptimize object that is the superset of the two values given as argument.
 */
constexpr inline ParameterToOptimize operator|(ParameterToOptimize p1, ParameterToOptimize p2) {
  using utype = std::underlying_type<ParameterToOptimize>::type;
  return static_cast<ParameterToOptimize>(static_cast<utype>(p1) | static_cast<utype>(p2));
}

class UpdateFunctionManager : public Utils::UpdateFunctionManagerBase {
 public:
  int hessianCounter;
  UpdateFunctionManager(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings,
                        ParameterToOptimize typeOfParameter);
  /**
   * @brief Update the errors vector for the least squares optimization for a given parameter
   */
  void updateErrors(const Eigen::VectorXd& parameters, Eigen::VectorXd& errors) override;
  /**
   * @brief This function returns the number of data points present in the least squares optimization,
   *        i.e., the number of Hessian elements for the parameter to be optimized and optionally
   *        2 data points to implement constraints on the parameter.
   *
   */
  int getNumberOfDataPoints(const Eigen::VectorXd& parameters) const override;
  /**
   * @brief Function that converts the variables of the optimization to the parameters object
   */
  MolecularMechanics::SfamParameters vectorToParameters(const Eigen::VectorXd& parameters);
  /**
   * @brief Function that converts the parameters object to the variables of the optimization
   */
  Eigen::VectorXd parametersToVector(const MolecularMechanics::SfamParameters& parameters);
  /**
   * @brief Set the initial variables of the optimization (variables are MM parameters).
   */
  void setInitialParameters(Eigen::VectorXd initialParameters);

 private:
  // The data used within all MM parametrization classes
  ParametrizationData& data_;
  // The settings
  std::shared_ptr<Utils::Settings> settings_;
  // Returns a list of atom indices for which the Hessian contributions need to be calculated.
  std::vector<int> getAtomsToConsiderForHessian() const;
  // The type of parameter that is going to be optimized here
  ParameterToOptimize typeOfParameter_;
  // The MM calculator
  std::unique_ptr<MolecularMechanics::SfamMolecularMechanicsCalculator> mmCalculator_;
  // Determines whether the given parameter type is a subset of the member typeOfParameter_
  bool typeOfParameterShouldBeOptimized(const ParameterToOptimize& typeToCheck) const {
    auto combined = typeToCheck | typeOfParameter_;
    return combined == typeOfParameter_;
  }
  // The initial parameters as an Eigen::VectorXd, this is saved so that the constraints can be set
  Eigen::VectorXd initialParameters_;
  // Whether to add constraints in the parameter optimization
  bool addConstraints_;
  // Holds the minimum and maximum values of the parameters as a scaling factor of the initial parameter.
  Eigen::VectorXd minimumParameters_;
  Eigen::VectorXd maximumParameters_;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_UPDATEFUNCTIONMANAGER_H
