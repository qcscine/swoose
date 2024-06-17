/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_DISPERSIONTERM_H
#define MOLECULARMECHANICS_DISPERSIONTERM_H

#include "Dispersion.h"
#include "InteractionTermBase.h"
#include <Utils/Dftd3/Dftd3.h>
#include <memory>

namespace Scine {

namespace Utils {
class AtomCollection;
class FullSecondDerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {
/**
 * @class DispersionTerm DispersionTerm.h
 * @brief Class evaluating dispersion (D3-BJ) interaction between two atoms.
 */
class DispersionTerm : public InteractionTermBase {
 public:
  using AtomIndex = int;

  /**
   * @brief Constructor.
   * @param firstAtom Index of first atom.
   * @param secondAtom Index of second atom.
   * @param dispersion An instance of the Dispersion class.
   */
  DispersionTerm(AtomIndex firstAtom, AtomIndex secondAtom, const Dispersion& dispersion, std::shared_ptr<double> cutoffRadius);
  /** @brief Destructor. */
  ~DispersionTerm();

  /**
   * @brief Evaluates energy contribution and adds the derivatives.
   */
  double evaluateDispersionTerm(const std::vector<Utils::Dftd3::Dftd3Atom>& structureOfDftd3Atoms,
                                Utils::FullSecondDerivativeCollection& derivatives,
                                std::shared_ptr<Utils::Dftd3::Dftd3> d3, Eigen::MatrixXd& R0) const;
  /** @brief Getter for index of first atom. */
  int getFirstAtom() const;
  /** @brief Getter for index of second atom. */
  int getSecondAtom() const;
  /** @brief Getter for the underlying instance of the Dispersion class. */
  Dispersion getDispersion() const;

 private:
  AtomIndex firstAtom_, secondAtom_;
  Dispersion dispersion_;
  // Pointer so that it is updated whenever the settings are updated.
  // Reference does not work because then the class loses its assignment operator.
  std::shared_ptr<double> cutoffRadius_;
};

} // namespace MolecularMechanics
} // namespace Scine
#endif // MOLECULARMECHANICS_DISPERSIONTERM_H
