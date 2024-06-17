/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_ANGLETERM_H
#define MOLECULARMECHANICS_ANGLETERM_H

#include "../Topology/AngleType.h"
#include "Angle.h"
#include "InteractionTermBase.h"
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {

namespace Utils {
class AtomicSecondDerivativeCollection;
} // namespace Utils

namespace MolecularMechanics {
/**
 * @class AngleTerm AngleTerm.h
 * @brief Class evaluating angle interaction for three given atoms.
 */
class AngleTerm : public InteractionTermBase {
 public:
  using AtomIndex = int;

  /** @brief Constructor from three atom types and instances of Angle and AngleType classes. */
  AngleTerm(AtomIndex firstAtom, AtomIndex secondAtom, AtomIndex thirdAtom, const Angle& angle, const AngleType& typeOfAngle);
  /** @brief Destructor. */
  ~AngleTerm();

  /**
   * @brief Evaluates energy contribution and adds the derivatives.
   */
  double evaluateAngleTerm(const Utils::PositionCollection& positions, Utils::AtomicSecondDerivativeCollection& derivatives) const;
  /**
   * @brief Getter for angle type.
   */
  AngleType getTypeOfAngle() const;

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

 private:
  friend class HydrogenBondTerm;
  static Utils::AutomaticDifferentiation::Second3D threeDimDer(const Utils::AutomaticDifferentiation::Second1D& energy,
                                                               const Utils::AutomaticDifferentiation::Second3D& alpha);
  static void calculateDerivativesWithNormalFormula(const Utils::AutomaticDifferentiation::Second1D& angleDerivative,
                                                    const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                                                    Utils::AutomaticDifferentiation::Second3D& angleContributionAtom1,
                                                    Utils::AutomaticDifferentiation::Second3D& angleContributionAtom2,
                                                    Utils::AutomaticDifferentiation::Second3D& angleContributionAtom3);
  static void calculateDerivativesForCriticalAngles(const Utils::AutomaticDifferentiation::Second1D& angleDerivative,
                                                    const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                                                    Utils::AutomaticDifferentiation::Second3D& angleContributionAtom1,
                                                    Utils::AutomaticDifferentiation::Second3D& angleContributionAtom2,
                                                    Utils::AutomaticDifferentiation::Second3D& angleContributionAtom3);

  AtomIndex firstAtom_, secondAtom_, thirdAtom_;
  Angle angle_;
  AngleType typeOfAngle_;
  static constexpr double singularityCriterion_ = 0.00001;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_ANGLETERM_H
