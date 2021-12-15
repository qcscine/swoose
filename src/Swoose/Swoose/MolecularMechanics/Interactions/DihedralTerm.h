/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_DIHEDRALTERM_H
#define MOLECULARMECHANICS_DIHEDRALTERM_H

#include "../Topology/DihedralType.h"
#include "Dihedral.h"
#include "InteractionTermBase.h"
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>
#include <Utils/Typenames.h>

namespace Scine {

namespace Utils {
class AtomicSecondDerivativeCollection;
class Tensor33;
} // namespace Utils

namespace MolecularMechanics {
/**
 * @class DihedralTerm DihedralTerm.h
 * @brief Class evaluating dihedral interaction for four given atoms.
 */
class DihedralTerm : public InteractionTermBase {
 public:
  using AtomIndex = int;

  /** @brief Constructor from four atom in indices and instances of Dihedral and DihedralType classes. */
  DihedralTerm(AtomIndex firstAtom, AtomIndex secondAtom, AtomIndex thirdAtom, AtomIndex fourthAtom,
               const Dihedral& dihedral, const DihedralType& typeOfDihedral);
  /** @brief Destructor. */
  ~DihedralTerm();

  /**
   * @brief Evaluates energy contribution and adds the derivatives.
   */
  double evaluateDihedralTerm(const Utils::PositionCollection& positions,
                              Utils::AtomicSecondDerivativeCollection& derivatives) const;

  /**
   * @brief Calculate and return the dihedral angle theta from the three essential vectors A, B and G describing it.
   */
  static double getTheta(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& G);
  /**
   * @brief Getter for the corresponding instance of the DihedralType class.
   */
  DihedralType getTypeOfDihedral() const;

  /**
   * @brief Getter for first atom.
   */
  int getFirstAtom() const;

  /**
   * @brief Getter for second atom.
   */
  int getSecondAtom() const;

  /**
   * @brief Getter for third atom.
   */
  int getThirdAtom() const;

  /**
   * @brief Getter for fourth atom.
   */
  int getFourthAtom() const;

 private:
  Utils::AutomaticDifferentiation::Second3D threeDimDer(const Utils::AutomaticDifferentiation::Second1D& energy,
                                                        const Utils::AutomaticDifferentiation::Second3D& alpha) const;
  void setSecondDerivative(Utils::AutomaticDifferentiation::Second3D& h, const Utils::Tensor33& tensor, double derFirst,
                           double derSecond) const;

  AtomIndex firstAtom_, secondAtom_, thirdAtom_, fourthAtom_;
  Dihedral dihedral_;
  DihedralType typeOfDihedral_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_DIHEDRALTERM_H