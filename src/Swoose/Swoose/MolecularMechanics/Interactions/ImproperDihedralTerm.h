/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMIMPROPERDIHEDRALTERM_H
#define MMIMPROPERDIHEDRALTERM_H

#include "../Topology/ImproperDihedralType.h"
#include "ImproperDihedral.h"
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
 * @class ImproperDihedralTerm ImproperDihedralTerm.h
 * @brief Class evaluating improper dihedral interaction for four given atoms.
 *
 *        This class is very similar to the DihedralTerm class.
 *        Here, instead of the third atom, we have the central atom.
 *        For future development: Combine classes in a smart way.
 */
class ImproperDihedralTerm : public InteractionTermBase {
 public:
  using AtomIndex = int;

  /**
   * @brief Constructor from four atoms and instances of an ImproperDihedral and an ImproperDihedralType.
   */
  ImproperDihedralTerm(AtomIndex firstAtom, AtomIndex secondAtom, AtomIndex centralAtom, AtomIndex fourthAtom,
                       const ImproperDihedral& improperDihedral, const ImproperDihedralType& typeOfImproperDihedral);
  /**
   * @brief Destructor.
   */
  ~ImproperDihedralTerm();

  /**
   * @brief Evaluates energy contribution and adds the derivatives.
   */
  double evaluateImproperDihedralTerm(const Utils::PositionCollection& positions,
                                      Utils::AtomicSecondDerivativeCollection& derivatives) const;

  /**
   * @brief Calculates and returns the improper dihedral angle theta from the three essential vectors A, B and G.
   */
  static double getTheta(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& G);

  /**
   * @brief Getter for the improper dihedral type.
   */
  ImproperDihedralType getTypeOfImproperDihedral() const;

  /**
   * @brief Getter for first atom.
   */
  int getFirstAtom() const;

  /**
   * @brief Getter for second atom.
   */
  int getSecondAtom() const;

  /**
   * @brief Getter for central atom.
   */
  int getCentralAtom() const;

  /**
   * @brief Getter for fourth atom.
   */
  int getFourthAtom() const;

 private:
  Utils::AutomaticDifferentiation::Second3D threeDimDer(const Utils::AutomaticDifferentiation::Second1D& energy,
                                                        const Utils::AutomaticDifferentiation::Second3D& alpha) const;
  void setSecondDerivative(Utils::AutomaticDifferentiation::Second3D& h, const Utils::Tensor33& tensor, double derFirst,
                           double derSecond) const;

  AtomIndex firstAtom_, secondAtom_, centralAtom_, fourthAtom_;
  ImproperDihedral improperDihedral_;
  ImproperDihedralType typeOfImproperDihedral_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MMIMPROPERDIHEDRALTERM_H
