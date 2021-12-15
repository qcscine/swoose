/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ParameterOptimizer.h"
#include "MMParametrizationSettings.h"
#include "OptimizationSetup.h"
#include "ParametrizationData.h"
#include "ParametrizationUtils/UpdateFunctionManager.h"
#include <Core/Log.h>
#include <Utils/Optimizer/LeastSquares/LevenbergMarquardt.h>

namespace {
constexpr const char* separator = "----------------------------------\n";
} // namespace

namespace Scine {
namespace MMParametrization {

ParameterOptimizer::ParameterOptimizer(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings, Core::Log& log)
  : data_(data), settings_(settings), log_(log) {
}

void ParameterOptimizer::optimizeParameters() {
  try {
    optimizeParametersImpl();
  }
  catch (std::exception& e) {
    throw std::runtime_error("The MM parameter optimization failed. Error: " + std::string(e.what()));
  }
}

void ParameterOptimizer::optimizeParametersImpl() {
  log_.output << separator;
  log_.output << "Starting to fit half barrier heights for dihedral potentials." << Core::Log::endl;

  // Fit dihedrals:
  optimizeDihedralParametersImpl();

  log_.output << separator;
  log_.output << "Starting to fit force constants for angle bending potentials." << Core::Log::endl;

  // Fit angles:
  optimizeAngleParametersImpl();

  log_.output << separator;
  log_.output << "Starting to fit force constants for bond stretching potentials." << Core::Log::endl;

  // Fit bonds:
  optimizeBondParametersImpl();

  // Fit improper dihedral force constants only if required
  if (settings_->getBool(SwooseUtilities::SettingsNames::optimizeImproperDihedralForceConstants)) {
    log_.output << separator;
    log_.output << "Starting to fit force constants for improper dihedral potentials." << Core::Log::endl;
    optimizeImproperDihedralParametersImpl();
  }

  log_.output << "MM parameter optimization done." << Core::Log::endl;
}

void ParameterOptimizer::optimizeDihedralParametersImpl() {
  auto dihedrals = data_.parameters.getDihedrals();
  for (const auto& dihedralParam : dihedrals) {
    if (dihedralParam.second.getHalfBarrierHeight() != OptimizationSetup::getInitialDihedralHalfBarrierHeight())
      continue;

    MolecularMechanics::DihedralType dihedralType1 = dihedralParam.first;
    log_.output << "Parameter description: " << dihedralType1.a1 << " " << dihedralType1.a2 << " " << dihedralType1.a3
                << " " << dihedralType1.a4 << Core::Log::nl;
    log_.output << "Value to optimize: " << dihedralParam.second.getHalfBarrierHeight() << " kcal/mol" << Core::Log::endl;
    data_.vectorOfAtomTypesForParameter = {dihedralType1.a1, dihedralType1.a2, dihedralType1.a3, dihedralType1.a4};

    for (const auto& dihedral : data_.topology.getDihedralContainer()) {
      auto atom1 = data_.atomTypes.getAtomType(dihedral.atom1);
      auto atom2 = data_.atomTypes.getAtomType(dihedral.atom2);
      auto atom3 = data_.atomTypes.getAtomType(dihedral.atom3);
      auto atom4 = data_.atomTypes.getAtomType(dihedral.atom4);
      MolecularMechanics::DihedralType dihedralType2(atom1, atom2, atom3, atom4);
      if (dihedralType2 == dihedralType1)
        data_.vectorOfHessianSubmatrixIndices.emplace_back(std::make_pair(dihedral.atom1, dihedral.atom4));
    }

    UpdateFunctionManager updateDihedral(data_, settings_, ParameterToOptimize::Dihedral);
    auto parametersAsVector = updateDihedral.parametersToVector(data_.parameters);
    updateDihedral.setInitialParameters(parametersAsVector);
    Utils::LevenbergMarquardt optDihedral;
    optDihedral.calculateCovarianceMatrix = false;
    optDihedral.maxFuncEval = maximumNumberOfFunctionEvaluationsInOptimization_;
    optDihedral.optimize(parametersAsVector, updateDihedral);
    log_.output << "This parameter has been optimized to be: " << parametersAsVector(0) << " kcal/mol" << Core::Log::nl;
    data_.parameters = updateDihedral.vectorToParameters(parametersAsVector);
    data_.vectorOfHessianSubmatrixIndices.clear(); // Empty the vector again for the next step
    log_.debug << "Hessian calculations: " << updateDihedral.hessianCounter << Core::Log::nl;
    log_.output << Core::Log::endl;
  }
}

void ParameterOptimizer::optimizeAngleParametersImpl() {
  auto angles = data_.parameters.getAngles();
  for (const auto& angleParam : angles) {
    if (angleParam.second.getForceConstant() != OptimizationSetup::getInitialAngleForceConstant())
      continue;

    MolecularMechanics::AngleType angleType1 = angleParam.first;
    log_.output << "Parameter description: " << angleType1.a1 << " " << angleType1.a2 << " " << angleType1.a3
                << Core::Log::nl;
    data_.vectorOfAtomTypesForParameter = {angleType1.a1, angleType1.a2, angleType1.a3};
    log_.output << "Value to optimize: " << angleParam.second.getForceConstant() << " kcal/(mol*rad^2)" << Core::Log::endl;

    for (const auto& angle : data_.topology.getAngleContainer()) {
      auto atom1 = data_.atomTypes.getAtomType(angle.atom1);
      auto atom2 = data_.atomTypes.getAtomType(angle.atom2);
      auto atom3 = data_.atomTypes.getAtomType(angle.atom3);
      MolecularMechanics::AngleType angleType2(atom1, atom2, atom3);
      if (angleType2 == angleType1)
        data_.vectorOfHessianSubmatrixIndices.emplace_back(std::make_pair(angle.atom1, angle.atom3));
    }

    UpdateFunctionManager updateAngle(data_, settings_, ParameterToOptimize::Angle);
    auto parametersAsVector = updateAngle.parametersToVector(data_.parameters);
    updateAngle.setInitialParameters(parametersAsVector);
    Utils::LevenbergMarquardt optAngle;
    optAngle.calculateCovarianceMatrix = false;
    optAngle.maxFuncEval = maximumNumberOfFunctionEvaluationsInOptimization_;
    optAngle.optimize(parametersAsVector, updateAngle);
    log_.output << "This parameter has been optimized to be: " << parametersAsVector(0) << " kcal/(mol*rad^2)"
                << Core::Log::nl;
    data_.parameters = updateAngle.vectorToParameters(parametersAsVector);
    data_.vectorOfHessianSubmatrixIndices.clear(); // Empty the vector again for the next step
    log_.debug << "Hessian calculations: " << updateAngle.hessianCounter << Core::Log::nl;
    log_.output << Core::Log::endl;
  }
}

void ParameterOptimizer::optimizeBondParametersImpl() {
  auto bonds = data_.parameters.getBonds();
  for (const auto& bondParam : bonds) {
    if (bondParam.second.getForceConstant() != OptimizationSetup::getInitialBondForceConstant())
      continue;

    MolecularMechanics::BondType bondType1 = bondParam.first;
    log_.output << "Parameter description: " << bondType1.a1 << " " << bondType1.a2 << Core::Log::nl;
    data_.vectorOfAtomTypesForParameter = {bondType1.a1, bondType1.a2};
    log_.output << "Value to optimize: " << bondParam.second.getForceConstant() << " kcal/(mol*A^2)" << Core::Log::endl;

    for (const auto& bond : data_.topology.getBondContainer()) {
      auto atom1 = data_.atomTypes.getAtomType(bond.atom1);
      auto atom2 = data_.atomTypes.getAtomType(bond.atom2);
      MolecularMechanics::BondType bondType2(atom1, atom2);
      if (bondType2 == bondType1) {
        data_.vectorOfHessianSubmatrixIndices.emplace_back(std::make_pair(bond.atom1, bond.atom2));
      }
    }

    UpdateFunctionManager updateBond(data_, settings_, ParameterToOptimize::Bond);
    auto parametersAsVector = updateBond.parametersToVector(data_.parameters);
    updateBond.setInitialParameters(parametersAsVector);
    Utils::LevenbergMarquardt optBond;
    optBond.calculateCovarianceMatrix = false;
    optBond.maxFuncEval = maximumNumberOfFunctionEvaluationsInOptimization_;
    optBond.optimize(parametersAsVector, updateBond);
    log_.output << "This parameter has been optimized to be: " << parametersAsVector(0) << " kcal/(mol*A^2)" << Core::Log::nl;
    data_.parameters = updateBond.vectorToParameters(parametersAsVector);
    data_.vectorOfHessianSubmatrixIndices.clear(); // Empty the vector again for the next step
    log_.debug << "Hessian calculations: " << updateBond.hessianCounter << Core::Log::nl;
    log_.output << Core::Log::endl;
  }
}

void ParameterOptimizer::optimizeImproperDihedralParametersImpl() {
  auto improperDihedrals = data_.parameters.getImproperDihedrals();
  for (const auto& improperDihedralParam : improperDihedrals) {
    // Check that it is not a planar group
    if (improperDihedralParam.second.getForceConstant() == OptimizationSetup::getImproperDihedralForceConstantForPlanarGroups())
      continue;

    // Check whether a value already exists for this parameter
    if (improperDihedralParam.second.getForceConstant() !=
        OptimizationSetup::getInitialImproperDihedralForceConstantForNonPlanarGroups())
      continue;

    MolecularMechanics::ImproperDihedralType improperDihedralType1 = improperDihedralParam.first;
    log_.output << "Parameter description: " << improperDihedralType1.ac << " " << improperDihedralType1.a2 << " "
                << improperDihedralType1.a3 << " " << improperDihedralType1.a4 << Core::Log::nl;
    log_.output << "Value to optimize: " << improperDihedralParam.second.getForceConstant() << " kcal/(mol*rad^2)"
                << Core::Log::endl;
    data_.vectorOfAtomTypesForParameter = {improperDihedralType1.ac, improperDihedralType1.a2, improperDihedralType1.a3,
                                           improperDihedralType1.a4};

    for (const auto& imp : data_.topology.getImproperDihedralContainer()) {
      auto atomC = data_.atomTypes.getAtomType(imp.centralAtom);
      auto atom2 = data_.atomTypes.getAtomType(imp.atom2);
      auto atom3 = data_.atomTypes.getAtomType(imp.atom3);
      auto atom4 = data_.atomTypes.getAtomType(imp.atom4);
      MolecularMechanics::ImproperDihedralType improperDihedralType2(atomC, atom2, atom3, atom4);
      if (improperDihedralType2 == improperDihedralType1) {
        data_.vectorOfHessianSubmatrixIndices.emplace_back(std::make_pair(imp.centralAtom, imp.centralAtom));
        data_.vectorOfHessianSubmatrixIndices.emplace_back(std::make_pair(imp.centralAtom, imp.atom2));
        data_.vectorOfHessianSubmatrixIndices.emplace_back(std::make_pair(imp.centralAtom, imp.atom3));
        data_.vectorOfHessianSubmatrixIndices.emplace_back(std::make_pair(imp.centralAtom, imp.atom4));
      }
    }

    UpdateFunctionManager updateImproperDihedral(data_, settings_, ParameterToOptimize::ImproperDihedral);
    auto parametersAsVector = updateImproperDihedral.parametersToVector(data_.parameters);
    updateImproperDihedral.setInitialParameters(parametersAsVector);
    Utils::LevenbergMarquardt optImproperDihedral;
    optImproperDihedral.calculateCovarianceMatrix = false;
    optImproperDihedral.maxFuncEval = maximumNumberOfFunctionEvaluationsInOptimization_;
    optImproperDihedral.optimize(parametersAsVector, updateImproperDihedral);
    log_.output << "This parameter has been optimized to be: " << parametersAsVector(0) << " kcal/(mol*rad^2)"
                << Core::Log::nl;
    data_.parameters = updateImproperDihedral.vectorToParameters(parametersAsVector);
    data_.vectorOfHessianSubmatrixIndices.clear(); // Empty the vector again for the next step
    log_.debug << "Hessian calculations: " << updateImproperDihedral.hessianCounter << Core::Log::nl;
    log_.output << Core::Log::endl;
  }
}

} // namespace MMParametrization
} // namespace Scine