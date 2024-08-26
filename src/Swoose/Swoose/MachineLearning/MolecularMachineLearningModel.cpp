/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MolecularMachineLearningModel.h"
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Math/MachineLearning/ChemicalRepresentations/AtomicForcesManager.h>
#include <Utils/Math/MachineLearning/ChemicalRepresentations/CoulombMatrix.h>
#include <Utils/Math/MachineLearning/CrossValidation.h>
#include <Utils/MolecularTrajectory.h>
#include <numeric>

namespace Scine {
namespace Swoose {
namespace MachineLearning {

void MolecularMachineLearningModel::setReferenceData(const Utils::MolecularTrajectory& structures,
                                                     const std::vector<ForcesCollection>& refForces) {
  if (refForces.size() != 0 && int(refForces.size()) != structures.size())
    throw std::runtime_error(
        "The number of structures is not equal to the number of individual collections of forces.");

  if (structures.empty())
    throw std::runtime_error("The molecular trajectory is empty.");

  const auto& elements = structures.getElementTypes();
  auto nAtoms = structures.molecularSize();

  Utils::AtomCollection originalStructure(elements, structures[0]);
  Utils::MachineLearning::AtomicForcesManager forcesManager(originalStructure);

  coulombMatrixFeatures_.clear();
  coulombMatrixFeatures_.reserve(structures.size());

  forcesFeatures_.clear();
  forcesFeatures_.reserve(structures.size());

  refEnergies_ = structures.getEnergies();
  refEnergies_.reserve(structures.size());

  refForces_.clear();
  refForces_.reserve(structures.size());

  forcePredictors_.resize(nAtoms);
  setDefaultHyperparameters(); // After resize of force predictors

  for (int structureIndex = 0; structureIndex < structures.size(); ++structureIndex) {
    // Handle energies
    if (!refEnergies_.empty()) {
      Utils::MachineLearning::CoulombMatrix cm({elements, structures[structureIndex]});
      coulombMatrixFeatures_.emplace_back(cm.getFeatures());
    }

    // Handle forces
    if (!refForces.empty()) {
      forcesManager.modifyPositions(structures[structureIndex]);
      std::vector<Eigen::VectorXd> features;
      ForcesCollection forces(nAtoms, 3);
      for (int i = 0; i < nAtoms; ++i) {
        features.emplace_back(forcesManager.calculateFeatures(i));
        Eigen::RowVector3d force = refForces[structureIndex].row(i);
        forces.row(i) = forcesManager.toInternalRepresentation(force, i);
      }
      refForces_.push_back(forces);
      forcesFeatures_.push_back(features);
    }
  }
}

void MolecularMachineLearningModel::trainEnergyModel() {
  if (coulombMatrixFeatures_.empty() || refEnergies_.empty())
    throw std::runtime_error("Cannot train energy model, data is missing.");

  auto features = getEnergyFeatures();
  auto targets = getEnergyTargets();

  // Train model
  energyPredictor_.trainModel(features, targets);
}

void MolecularMachineLearningModel::trainForcesModel() {
  if (forcesFeatures_.empty() || refForces_.empty())
    throw std::runtime_error("Cannot train forces model, data is missing.");

  auto nAtoms = refForces_[0].rows();

#pragma omp parallel for
  for (int i = 0; i < nAtoms; ++i) {
    trainSingleForceModel(i);
  }
}

void MolecularMachineLearningModel::trainSingleForceModel(int atomIndex) {
  // Get features and targets for the atomic force of the atom with the index 'atomIndex'
  auto features = getSingleForceFeatures(atomIndex);
  auto targets = getSingleForceTargets(atomIndex);

  // Train model
  forcePredictors_[atomIndex].trainModel(features, targets);
}

double MolecularMachineLearningModel::predictEnergy(const Utils::AtomCollection& structure) {
  Utils::MachineLearning::CoulombMatrix cm(structure);
  auto prediction = energyPredictor_.predict(cm.getFeatures());
  return prediction(0);
}

ForcesCollection MolecularMachineLearningModel::predictForces(const Utils::AtomCollection& structure) {
  assert(refForces_.at(0).rows() == structure.size());
  Utils::MachineLearning::AtomicForcesManager forcesManager(structure);
  auto nAtoms = refForces_.at(0).rows();
  ForcesCollection forces(nAtoms, 3);
  for (int i = 0; i < refForces_.at(0).rows(); ++i) {
    Eigen::RowVector3d prediction = forcePredictors_[i].predict(forcesManager.calculateFeatures(i));
    forces.row(i) = forcesManager.toGlobalRepresentation(prediction, i);
  }
  return forces;
}

std::pair<double, double> MolecularMachineLearningModel::evaluateEnergyModel(int k) {
  if (coulombMatrixFeatures_.empty() || refEnergies_.empty())
    throw std::runtime_error("Cannot validate the energy model, data is missing.");

  auto features = getEnergyFeatures();
  auto targets = getEnergyTargets();

  Utils::MachineLearning::CrossValidation cv(energyPredictor_, k);
  return cv.evaluateRegressionModel(features, targets);
}

std::pair<double, double> MolecularMachineLearningModel::evaluateForcesModel(int k, bool pooledVariance) {
  if (forcesFeatures_.empty() || refForces_.empty())
    throw std::runtime_error("Cannot validate the atomic forces model, data is missing.");

  auto nAtoms = refForces_[0].rows();
  std::vector<double> variances(nAtoms);
  std::vector<double> errors(nAtoms);

#pragma omp parallel for
  for (int i = 0; i < nAtoms; ++i) {
    auto result = evaluateSingleForceModel(i, k);
    errors[i] = result.first;
    variances[i] = result.second * result.second;
  }

  // Calculate average MAE
  auto combinedMeanAbsoluteError = std::accumulate(errors.begin(), errors.end(), 0.0) / errors.size();

  if (pooledVariance)
    return std::make_pair(combinedMeanAbsoluteError,
                          sqrt(std::accumulate(variances.begin(), variances.end(), 0.0) / variances.size()));

  // Calculate true combined variance:
  auto total = 0.0;
  for (int i = 0; i < int(variances.size()); ++i) {
    total += variances[i] + (errors[i] - combinedMeanAbsoluteError) * (errors[i] - combinedMeanAbsoluteError);
  }

  return std::make_pair(combinedMeanAbsoluteError, sqrt(total / variances.size()));
}

std::pair<double, double> MolecularMachineLearningModel::evaluateSingleForceModel(int atomIndex, int k) {
  auto features = getSingleForceFeatures(atomIndex);
  auto targets = getSingleForceTargets(atomIndex);
  Utils::MachineLearning::CrossValidation cv(forcePredictors_[atomIndex], k);
  return cv.evaluateRegressionModel(features, targets);
}

Eigen::MatrixXd MolecularMachineLearningModel::getEnergyFeatures() {
  Eigen::MatrixXd features(coulombMatrixFeatures_.size(), coulombMatrixFeatures_[0].size());
  for (int i = 0; i < features.rows(); ++i) {
    features.row(i) = coulombMatrixFeatures_[i];
  }
  return features;
}

Eigen::VectorXd MolecularMachineLearningModel::getEnergyTargets() {
  return Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(refEnergies_.data(), refEnergies_.size());
}

Eigen::MatrixXd MolecularMachineLearningModel::getSingleForceFeatures(int atomIndex) {
  int nDataPoints = forcesFeatures_.size();
  Eigen::MatrixXd features(nDataPoints, forcesFeatures_[0].at(atomIndex).size());
  for (int i = 0; i < nDataPoints; ++i) {
    features.row(i) = forcesFeatures_[i].at(atomIndex);
  }
  return features;
}

Eigen::MatrixXd MolecularMachineLearningModel::getSingleForceTargets(int atomIndex) {
  const size_t nDataPoints = forcesFeatures_.size();
  Eigen::MatrixXd targets(nDataPoints, 3);
  for (size_t j = 0; j < nDataPoints; ++j) {
    targets.row(j) = refForces_[j].row(atomIndex);
  }
  return targets;
}

Utils::MachineLearning::KernelRidgeRegression& MolecularMachineLearningModel::energyPredictor() {
  return energyPredictor_;
}

Utils::MachineLearning::KernelRidgeRegression& MolecularMachineLearningModel::forcePredictor(int atomIndex) {
  return forcePredictors_.at(atomIndex);
}

// TODO: The hyperparameters should be optimized again.
// TODO: The choice of the kernel should be depending on amount of available reference data.
void MolecularMachineLearningModel::setDefaultHyperparameters() {
  Eigen::VectorXd sigma(1);
  sigma[0] = 50.0;
  energyPredictor_.setKernel(Utils::MachineLearning::Kernels::laplacianKernel, sigma);
  for (auto& predictor : forcePredictors_) {
    predictor.setKernel(Utils::MachineLearning::Kernels::laplacianKernel, sigma);
  }
}

} // namespace MachineLearning
} // namespace Swoose
} // namespace Scine
