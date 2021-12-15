/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMIMPROPERDIHEDRAL_H
#define MMIMPROPERDIHEDRAL_H

#include <Utils/Math/AutomaticDifferentiation/Second1D.h>

namespace Scine {
namespace MolecularMechanics {
/**
 * @class ImproperDihedral ImproperDihedral.h
 * @brief Class treating an improper dihedral interaction, based solely on the angle (in rad), i.e. in 1 dimension.
 */
class ImproperDihedral {
 public:
  /**
   * @brief Constructor from a force constant and an equilibrium angle.
   */
  ImproperDihedral(double forceConstant, double equilibriumAngle);
  /**
   * @brief Constructor without arguments sets both parameters to zero and records that no parameters
   *        are available for this improper dihedral, hence throwing an exception if getInteraction() is called.
   *        This will not throw any error if the corresponding improper dihedral term
   *        is disabled (e.g., in a QM/MM calculation).
   */
  ImproperDihedral();
  /**
   * @brief Returns whether parameters are available for this improper dihedral.
   */
  bool hasParameters() const;
  /**
   * @brief Calculates the energy and its derivative for a given angle.
   */
  Utils::AutomaticDifferentiation::Second1D getInteraction(double angle) const;
  /**
   * @brief Getter for the threshold that decides whether a group is considered planar for improper dihedrals.
   */
  static double getPlanarGroupThreshold();

 private:
  double forceConstant_;
  double equilibriumAngle_;
  bool parametersAreAvailable_;
  static constexpr double planarGroupThreshold_ = 20.0;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MMIMPROPERDIHEDRAL_H