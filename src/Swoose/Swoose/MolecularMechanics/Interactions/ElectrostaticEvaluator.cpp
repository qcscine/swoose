/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Swoose/MolecularMechanics/Interactions/ElectrostaticEvaluator.h"
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Scine {
namespace MolecularMechanics {

ElectrostaticEvaluator::ElectrostaticEvaluator(const Utils::PositionCollection& positions, const std::vector<double>& atomicCharges)
  : InteractionExclusion(positions.rows()),
    ScaledInteractions(positions.rows()),
    positions_(positions),
    atomicCharges_(atomicCharges) {
}

double ElectrostaticEvaluator::evaluate(Utils::DerivativeCollection& derivatives) {
  const unsigned int nAtoms = this->positions_.rows();
  const Eigen::Matrix<bool, Eigen::Dynamic, 1> allAtoms = Eigen::Matrix<bool, Eigen::Dynamic, 1>::Constant(nAtoms, 1, true);
  const Eigen::SparseMatrix<bool>& exclusions = this->getExclusions();
  const Eigen::SparseMatrix<bool>& scaledTerms = this->getScaledInteractionPairs();
  const double& scalingFactor = this->getInteractionScalingFactor();

  // Avoid nested parallelization.
  omp_set_max_active_levels(1);
  const unsigned int nThreads = omp_get_max_threads();
  Eigen::VectorXd energySet = Eigen::VectorXd::Zero(nThreads);
  std::vector<Utils::DerivativeCollection> derivativeSet(nThreads,
                                                         Utils::DerivativeCollection(int(nAtoms), derivatives.getOrder()));
  for (auto& deriv : derivativeSet) {
    deriv.setZero();
  }
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iAtom = 0; iAtom < nAtoms; ++iAtom) {
    const unsigned int threadID = omp_get_thread_num();
    // All atoms for which the interaction is neither scaled nor excluded.
    // We must call the method pruned() here to ensure that any entry with false is removed from the sparse vector.
    Eigen::SparseVector<bool> notScaledOrExcluded =
        (allAtoms - (exclusions.col(iAtom) || scaledTerms.col(iAtom))).eval().pruned();
    // Evaluate non-scaled terms.
    energySet(threadID) += evaluateTermsForAtom(iAtom, derivativeSet[threadID], notScaledOrExcluded, 1.0);
    // Evaluate scaled terms.
    const Eigen::SparseVector<bool> scaledButExcluded = scaledTerms.col(iAtom) && exclusions.col(iAtom);
    energySet(threadID) += evaluateTermsForAtom(iAtom, derivativeSet[threadID],
                                                (scaledTerms.col(iAtom) - scaledButExcluded).pruned(), scalingFactor);
  } // for iAtom
  for (const auto& deriv : derivativeSet) {
    derivatives += deriv;
  }
  return energySet.sum();
}

double ElectrostaticEvaluator::evaluateTermsForAtom(unsigned int atomIndex, Utils::DerivativeCollection& derivatives,
                                                    const Eigen::SparseVector<bool>& otherAtoms, double scaling) {
  double energyIncrement = 0.0;
  const auto positionsI = this->positions_.row(atomIndex);
  const double& chargeI = this->atomicCharges_[atomIndex];
  const double totalScalingFactor = scalingFactor_ * scalingFactor_ * scaling;
  for (Eigen::SparseVector<bool>::InnerIterator it(otherAtoms); it; ++it) {
    if (it.row() >= atomIndex) {
      break;
    }
    const unsigned int atomIndexII = it.row();
    const auto positionDifference = this->positions_.row(atomIndexII) - positionsI;
    const double distance = positionDifference.norm();
    if (distance > *this->cutOffRadius_) {
      continue;
    }
    const double& chargeII = this->atomicCharges_[atomIndexII];
    const auto contribution =
        totalScalingFactor /
        Utils::AutomaticDifferentiation::variableWithUnitDerivative<Utils::DerivativeOrder::Two>(distance) *
        (chargeI * chargeII);
    const auto value =
        Utils::AutomaticDifferentiation::get3Dfrom1D<Utils::DerivativeOrder::Two>(contribution, positionDifference);
    derivatives.addDerivative(int(atomIndex), int(atomIndexII), value);
    energyIncrement += value.value();
  }
  return energyIncrement;
}

const std::vector<double>& ElectrostaticEvaluator::getAtomicCharges() {
  return atomicCharges_;
}

void ElectrostaticEvaluator::setScalingFactor(const double& scalingFactor) {
  scalingFactor_ = scalingFactor;
}
void ElectrostaticEvaluator::setCutOffRadius(std::shared_ptr<double> cutOffRadius) {
  cutOffRadius_ = cutOffRadius;
}

} // namespace MolecularMechanics
} // namespace Scine
