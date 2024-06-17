/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMMMHESSIANEVALUATOR_H
#define SWOOSE_QMMM_QMMMHESSIANEVALUATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/DataStructures/PartialHessian.h>
#include <memory>

namespace Scine {

namespace Utils {
class AtomCollection;
} // namespace Utils

namespace Qmmm {

class QmmmHessianEvaluator {
 public:
  /**
   * @brief Constructor.
   * @param qmCalculator
   * @param listOfQmAtoms A vector containing the indices of the QM atoms.
   */
  QmmmHessianEvaluator(std::shared_ptr<Core::Calculator> qmCalculator, const std::vector<int>& listOfQmAtoms);

  /**
   * @brief Calculates the Partial Hessian matrix from all the components given in the constructor.
   * @return The QM Hessian.
   */
  Utils::PartialHessian calculatePartialHessian();

 private:
  std::shared_ptr<Core::Calculator> qmCalculator_;
  const std::vector<int>& listOfQmAtoms_;
};

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMMMHESSIANEVALUATOR_H
