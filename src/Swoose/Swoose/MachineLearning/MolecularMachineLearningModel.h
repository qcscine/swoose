/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_MOLECULARMACHINELEARNINGMODEL_H
#define SWOOSE_MOLECULARMACHINELEARNINGMODEL_H

#include <Utils/Math/MachineLearning/Regression/KernelRidgeRegression.h>
#include <Eigen/Dense>
#include <vector>

namespace Scine {

// Forward declarations
namespace Utils {
class MolecularTrajectory;
class AtomCollection;
} // namespace Utils

namespace Swoose {
using ForcesCollection = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

namespace MachineLearning {

/**
 * @class MolecularMachineLearningModel MolecularMachineLearningModel.h
 * @brief The combination of the machine learning models for molecular energies and atomic forces.
 */
class MolecularMachineLearningModel {
 public:
  /**
   * @brief Default constructor.
   */
  MolecularMachineLearningModel() = default;
  /**
   * @brief Sets the reference data.
   * @param structures The molecular trajectory containing the structures along with their reference energies.
   * @param refForces A vector of the reference forces for each structure.
   */
  void setReferenceData(const Utils::MolecularTrajectory& structures, const std::vector<ForcesCollection>& refForces);
  /**
   * @brief Trains the energy model.
   */
  void trainEnergyModel();
  /**
   * @brief Trains the forces model.
   */
  void trainForcesModel();
  /**
   * @brief Predicts the energy for a given structure.
   */
  double predictEnergy(const Utils::AtomCollection& structure);
  /**
   * @brief Predicts the forces for a given structure.
   */
  ForcesCollection predictForces(const Utils::AtomCollection& structure);
  /**
   * @brief Validates the energy model via k-fold cross validation.
   * @param k The number of subsets k for the k-fold cross validation algorithm.
   * @return The mean absolute error and the standard deviation of the model evaluation as a std::pair.
   */
  std::pair<double, double> evaluateEnergyModel(int k);
  /**
   * @brief Validates the model for the atomic forces via k-fold cross validation.
   * @param k The number of subsets k for the k-fold cross validation algorithm.
   * @param pooledVariance Decides whether the standard deviation estimate for the forces
   *                       is based on the principle of pooled variance. The alternative is a true
   *                       combined variance that considers the different MAE means of the atomic forces.
   * @return The combined mean absolute error and the combined standard deviation
   *         of the model evaluation for the individual atomic forces as a std::pair.
   */
  std::pair<double, double> evaluateForcesModel(int k, bool pooledVariance = true);
  /**
   * @brief Accessor for the underlying energy model.
   */
  Utils::MachineLearning::KernelRidgeRegression& energyPredictor();
  /**
   * @brief Accessor for the underlying force model of the atom with index 'atomIndex'.
   */
  Utils::MachineLearning::KernelRidgeRegression& forcePredictor(int atomIndex);

 private:
  // Sets the default values of the model hyperparameters for all of the predictors.
  void setDefaultHyperparameters();
  // Trains a single forces model for a given atom index
  void trainSingleForceModel(int atomIndex);
  // Evaluates a single force model for a given atom index with k-fold cross-validation.
  std::pair<double, double> evaluateSingleForceModel(int atomIndex, int k);
  // Returns the energy feature matrix
  Eigen::MatrixXd getEnergyFeatures();
  // Returns the energy target vector
  Eigen::VectorXd getEnergyTargets();
  // Returns the feature matrix for the atomic force of the atom with index 'atomIndex'
  Eigen::MatrixXd getSingleForceFeatures(int atomIndex);
  // Returns the target matrix for the atomic force of the atom with index 'atomIndex'
  Eigen::MatrixXd getSingleForceTargets(int atomIndex);
  // The KRR model that learns and predicts the energy of the molecular system
  Utils::MachineLearning::KernelRidgeRegression energyPredictor_;
  // The KRR models that learn and predict the atomic forces of the molecular system.
  // It exists one model for each atom.
  std::vector<Utils::MachineLearning::KernelRidgeRegression> forcePredictors_;
  // The Coulomb matrix features for all structures (reference data)
  std::vector<Eigen::VectorXd> coulombMatrixFeatures_;
  // The features for learning atomic forces for all structures
  std::vector<std::vector<Eigen::VectorXd>> forcesFeatures_;
  // Reference energy targets
  std::vector<double> refEnergies_;
  // Reference atomic forces targets
  std::vector<ForcesCollection> refForces_;
};

} // namespace MachineLearning
} // namespace Swoose
} // namespace Scine

#endif // SWOOSE_MOLECULARMACHINELEARNINGMODEL_H
