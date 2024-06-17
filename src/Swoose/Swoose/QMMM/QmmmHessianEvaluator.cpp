/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmmmHessianEvaluator.h"
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry/AtomCollection.h>
#include <utility>

namespace Scine {
namespace Qmmm {

QmmmHessianEvaluator::QmmmHessianEvaluator(std::shared_ptr<Core::Calculator> qmCalculator, const std::vector<int>& listOfQmAtoms)
  : qmCalculator_(std::move(qmCalculator)), listOfQmAtoms_(listOfQmAtoms) {
}

Utils::PartialHessian QmmmHessianEvaluator::calculatePartialHessian() {
  Utils::PropertyList originalProperties = qmCalculator_->getRequiredProperties();
  qmCalculator_->setRequiredProperties(Utils::Property::Energy | Utils::Property::Hessian);
  auto result = qmCalculator_->calculate("");
  qmCalculator_->setRequiredProperties(originalProperties);

  Utils::HessianMatrix hessianMatrix = result.get<Utils::Property::Hessian>();
  const int numQmAtoms = listOfQmAtoms_.size();
  // this works because the link atoms are always at the end of the AtomCollection.
  Eigen::MatrixXd partialHessianMatrix = hessianMatrix.block(0, 0, 3 * numQmAtoms, 3 * numQmAtoms);
  return Utils::PartialHessian(partialHessianMatrix, listOfQmAtoms_);
}

} // namespace Qmmm
} // namespace Scine
