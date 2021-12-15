/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "UpdateFunctionManager.h"
#include "../MMParametrizationSettings.h"
#include "../ParametrizationData.h"
#include <Swoose/MolecularMechanics/SFAM/SfamCalculatorSettings.h>

namespace Scine {
namespace MMParametrization {

// This number controls the severity of the constraints
static constexpr double severityFactorForConstraints = 10000.0; // TODO: Adapt to taking away the scaling factor below.

UpdateFunctionManager::UpdateFunctionManager(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings,
                                             ParameterToOptimize typeOfParameter)
  : hessianCounter(0),
    data_(data),
    settings_(settings),
    typeOfParameter_(typeOfParameter),
    mmCalculator_(std::make_unique<MolecularMechanics::SfamMolecularMechanicsCalculator>()) {
  addConstraints_ = settings_->getBool(SwooseUtilities::SettingsNames::constrainMMParameters);
  auto silentLogger = Core::Log::silent();
  mmCalculator_->setLog(silentLogger);
  mmCalculator_->settings().modifyBool(SwooseUtilities::SettingsNames::applyCutoffDuringInitialization, true);
  mmCalculator_->settings().modifyBool(SwooseUtilities::SettingsNames::hydrogenBondCorrection, false);  // TODO
  mmCalculator_->settings().modifyDouble(SwooseUtilities::SettingsNames::nonCovalentCutoffRadius, 6.0); // TODO
  mmCalculator_->setListsOfNeighbors(data_.listsOfNeighbors);
  mmCalculator_->setParameters(data_.parameters);
  mmCalculator_->setStructure(data_.fullStructure);
  mmCalculator_->setRequiredProperties(Utils::Property::Hessian);
}

void UpdateFunctionManager::updateErrors(const Eigen::VectorXd& parameters, Eigen::VectorXd& errors) {
  auto numberOfParameters = static_cast<int>(parameters.size());
  minimumParameters_.resize(numberOfParameters);
  maximumParameters_.resize(numberOfParameters);

  // Get SfamParameters object from Eigen::VectorXd
  auto mmParameters = vectorToParameters(parameters);

  // Resize the error vector that will be updated in this function
  errors.resize(getNumberOfDataPoints(parameters));
  errors.setZero();

  // Counter variable for the errors vector
  int i = 0;

  // Prepare MM calculator
  mmCalculator_->setParameters(mmParameters); // for the cloned calculators in the Hessian calculator
  mmCalculator_->generatePotentialTerms(mmParameters, data_.topology, data_.atomTypes);

  mmCalculator_->setRequiredProperties(Utils::Property::Hessian); // TODO: Why is this needed, bug?

  mmCalculator_->setAtomsToConsiderForHessian(getAtomsToConsiderForHessian());
  mmCalculator_->calculate("");
  hessianCounter++;
  Utils::HessianMatrix mmHessianMatrix = mmCalculator_->results().take<Utils::Property::Hessian>();

  for (const auto& hessianBlockIndices : data_.vectorOfHessianSubmatrixIndices) {
    const auto& hessianSubmatrix = mmHessianMatrix.block<3, 3>(3 * hessianBlockIndices.first, 3 * hessianBlockIndices.second);
    auto referenceHessianSubmatrix =
        data_.fullHessian.block(3 * hessianBlockIndices.first, 3 * hessianBlockIndices.second, 3, 3);

    for (int k = 0; k < hessianSubmatrix.rows(); ++k) {
      for (int l = 0; l < hessianSubmatrix.cols(); ++l) {
        double hessianDifference = hessianSubmatrix(k, l) - referenceHessianSubmatrix.coeffRef(k, l);
        // TODO: REMOVE THE SCALING IN THE NEXT LINE EVENTUALLY, BUT FIRST CHECK WHETHER lambda HAS TO BE ADAPTED THEN
        hessianDifference *= Utils::Constants::kCalPerMol_per_hartree;
        errors(i) = hessianDifference;
        i++;
      }
    }
  }

  if (addConstraints_) {
    for (int j = 0; j < numberOfParameters; ++j) {
      double penaltyMin = 0.0;
      double penaltyMax = 0.0;
      if (minimumParameters_(j) < 1.0) {
        penaltyMin = std::min(0.0, parameters(j) - initialParameters_(j) * minimumParameters_(j));
      }
      if (maximumParameters_(j) > 1.0) {
        penaltyMax = std::max(0.0, parameters(j) - initialParameters_(j) * maximumParameters_(j));
      }
      errors(i) = severityFactorForConstraints * penaltyMin;
      i++;
      errors(i) = severityFactorForConstraints * penaltyMax;
      i++;
    }
  }
  assert(i == getNumberOfDataPoints(parameters));
}

// Conversion of the Eigen::VectorXd holding the parameters to optimize and the SfamParameters object
MolecularMechanics::SfamParameters UpdateFunctionManager::vectorToParameters(const Eigen::VectorXd& parameters) {
  MolecularMechanics::SfamParameters mmParameters = data_.parameters; // TODO: Is this copy really necessary here?

  // Counter for the parameters vector
  int i = 0;

  // For bonds
  if (typeOfParameterShouldBeOptimized(ParameterToOptimize::Bond)) {
    MolecularMechanics::BondType bondType(data_.vectorOfAtomTypesForParameter[0], data_.vectorOfAtomTypesForParameter[1]);
    for (auto& bond : mmParameters.getBonds()) {
      if (bond.first == bondType) {
        bond.second.setForceConstant(parameters(i));
        minimumParameters_(i) = 0.05; // maximum decrease with respect to initial parameter
        maximumParameters_(i) = 5.0;  // maximum increase with respect to initial parameter
        i++;
      }
    }
  }

  // For angles
  if (typeOfParameterShouldBeOptimized(ParameterToOptimize::Angle)) {
    MolecularMechanics::AngleType angleType(data_.vectorOfAtomTypesForParameter[0], data_.vectorOfAtomTypesForParameter[1],
                                            data_.vectorOfAtomTypesForParameter[2]);
    for (auto& angle : mmParameters.getAngles()) {
      if (angle.first == angleType) {
        angle.second.setForceConstant(parameters(i));
        minimumParameters_(i) = 0.05; // maximum decrease with respect to initial parameter
        maximumParameters_(i) = 5.0;  // maximum increase with respect to initial parameter
        i++;
      }
    }
  }

  // For dihedrals
  if (typeOfParameterShouldBeOptimized(ParameterToOptimize::Dihedral)) {
    MolecularMechanics::DihedralType dihedralType(
        data_.vectorOfAtomTypesForParameter[0], data_.vectorOfAtomTypesForParameter[1],
        data_.vectorOfAtomTypesForParameter[2], data_.vectorOfAtomTypesForParameter[3]);
    for (auto& dihedral : mmParameters.getDihedrals()) {
      if (dihedral.first == dihedralType) {
        dihedral.second.setHalfBarrierHeight(parameters(i));
        minimumParameters_(i) = 0.0; // decrease constraint to zero, i.e., no negative values allowed
        maximumParameters_(i) = 1.0; // no constraint
        i++;
      }
    }
  }

  // For improper dihedrals
  if (typeOfParameterShouldBeOptimized(ParameterToOptimize::ImproperDihedral)) {
    MolecularMechanics::ImproperDihedralType improperType(
        data_.vectorOfAtomTypesForParameter[0], data_.vectorOfAtomTypesForParameter[1],
        data_.vectorOfAtomTypesForParameter[2], data_.vectorOfAtomTypesForParameter[3]);
    for (auto& improperDihedral : mmParameters.getImproperDihedrals()) {
      if (improperDihedral.first == improperType) {
        improperDihedral.second.setForceConstant(parameters(i));
        minimumParameters_(i) = 1.0; // no constraint // TODO: Find the right bounds.
        maximumParameters_(i) = 1.0; // no constraint
        i++;
      }
    }
  }

  return mmParameters;
}

// Conversion of the parameter object to an Eigen::VectorXd of the parameters for the least squares optimization
Eigen::VectorXd UpdateFunctionManager::parametersToVector(const MolecularMechanics::SfamParameters& parameters) {
  std::vector<double> variables; // Start with std::vector such that we can use the push_back function.

  // For bonds
  if (typeOfParameterShouldBeOptimized(ParameterToOptimize::Bond)) {
    MolecularMechanics::BondType bondType(data_.vectorOfAtomTypesForParameter[0], data_.vectorOfAtomTypesForParameter[1]);
    for (const auto& bond : parameters.getBonds()) {
      if (bond.first == bondType)
        variables.push_back(bond.second.getForceConstant());
    }
  }

  // For angles
  if (typeOfParameterShouldBeOptimized(ParameterToOptimize::Angle)) {
    MolecularMechanics::AngleType angleType(data_.vectorOfAtomTypesForParameter[0], data_.vectorOfAtomTypesForParameter[1],
                                            data_.vectorOfAtomTypesForParameter[2]);
    for (const auto& angle : parameters.getAngles()) {
      if (angle.first == angleType)
        variables.push_back(angle.second.getForceConstant());
    }
  }

  // For dihedrals
  if (typeOfParameterShouldBeOptimized(ParameterToOptimize::Dihedral)) {
    MolecularMechanics::DihedralType dihedralType(
        data_.vectorOfAtomTypesForParameter[0], data_.vectorOfAtomTypesForParameter[1],
        data_.vectorOfAtomTypesForParameter[2], data_.vectorOfAtomTypesForParameter[3]);
    for (const auto& dihedral : parameters.getDihedrals()) {
      if (dihedral.first == dihedralType)
        variables.push_back(dihedral.second.getHalfBarrierHeight());
    }
  }

  // For improper dihedrals
  if (typeOfParameterShouldBeOptimized(ParameterToOptimize::ImproperDihedral)) {
    MolecularMechanics::ImproperDihedralType improperType(
        data_.vectorOfAtomTypesForParameter[0], data_.vectorOfAtomTypesForParameter[1],
        data_.vectorOfAtomTypesForParameter[2], data_.vectorOfAtomTypesForParameter[3]);
    for (const auto& improperDihedral : parameters.getImproperDihedrals()) {
      if (improperDihedral.first == improperType)
        variables.push_back(improperDihedral.second.getForceConstant());
    }
  }

  // Convert std::vector<double> to Eigen::VectorXd and return it
  Eigen::VectorXd eigenObject = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(variables.data(), variables.size());
  return eigenObject;
}

int UpdateFunctionManager::getNumberOfDataPoints(const Eigen::VectorXd& parameters) const {
  auto numberOfParameters = static_cast<int>(parameters.size());
  int dataPoints = 0; // initialize to zero
  auto numberOfHessianSubmatrices = static_cast<int>(data_.vectorOfHessianSubmatrixIndices.size());
  dataPoints += 9 * numberOfHessianSubmatrices; // 9 elements in each partial Hessian submatrix

  if (addConstraints_)
    dataPoints += 2 * numberOfParameters;

  return dataPoints;
}

void UpdateFunctionManager::setInitialParameters(Eigen::VectorXd initialParameters) {
  initialParameters_ = std::move(initialParameters);
}

std::vector<int> UpdateFunctionManager::getAtomsToConsiderForHessian() const {
  std::vector<int> result;
  for (const auto& hessianBlockIndices : data_.vectorOfHessianSubmatrixIndices) {
    if (std::find(result.begin(), result.end(), hessianBlockIndices.first) == result.end())
      result.push_back(hessianBlockIndices.first);
    if (std::find(result.begin(), result.end(), hessianBlockIndices.second) == result.end())
      result.push_back(hessianBlockIndices.second);
  }
  return result;
}

} // namespace MMParametrization
} // namespace Scine
