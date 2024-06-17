/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_ELECTROSTATICTERM_H
#define MOLECULARMECHANICS_ELECTROSTATICTERM_H

#include "Electrostatic.h"
#include "InteractionTermBase.h"
#include <Utils/Typenames.h>
#include <memory>
#include <vector>

namespace Scine {

namespace Utils {
class FullSecondDerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {
/**
 * @class ElectrostaticTerm ElectrostaticTerm.h
 * @brief Class evaluating electrostatic interaction between two atoms.
 */
class ElectrostaticTerm : public InteractionTermBase {
 public:
  using AtomIndex = int;

  /**
   * @brief Constructor from two atom indices and an instance of the Electrostatic class.
   */
  ElectrostaticTerm(AtomIndex firstAtom, AtomIndex secondAtom, const Electrostatic& electrostatic,
                    std::shared_ptr<double> cutoffRadius);

  /**
   * @brief Evaluates energy contribution and adds the derivatives.
   */
  double evaluateElectrostaticTerm(const Utils::PositionCollection& positions, Utils::FullSecondDerivativeCollection& derivatives,
                                   const std::vector<double>& atomicCharges, const double& scalingFactorForEachCharge) const;

  /** @brief Getter for index of first atom. */
  int getFirstAtom() const;
  /** @brief Getter for index of second atom. */
  int getSecondAtom() const;

 private:
  AtomIndex firstAtom_, secondAtom_;
  Electrostatic electrostatic_;
  // Pointer so that it is updated whenever the settings are updated.
  // Reference does not work because then the class loses its assignment operator.
  std::shared_ptr<double> cutoffRadius_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_ELECTROSTATICTERM_H
