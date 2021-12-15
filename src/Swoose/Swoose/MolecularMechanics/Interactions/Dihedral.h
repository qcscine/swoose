/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_DIHEDRAL_H
#define MOLECULARMECHANICS_DIHEDRAL_H

#include <Utils/Math/AutomaticDifferentiation/Second1D.h>

namespace Scine {
namespace MolecularMechanics {
/**
 * @class Dihedral Dihedral.h
 * @brief Class treating a dihedral interaction, based solely on the angle (in rad), i.e. in 1 dimension.
 */
class Dihedral {
 public:
  /** @brief Constructor from the three dihedral potential parameters. */
  Dihedral(double halfBarrierHeight, int periodicity, double phaseShift);
  /**
   * @brief Constructor without arguments sets all three parameters to zero and records that no parameters
   *        are available for this dihedral, hence throwing an exception if getInteraction() is called. This
   *        will not throw any error if the corresponding bonded term is disabled (e.g., in a QM/MM calculation).
   */
  Dihedral();
  /** @brief Calculates the energy with its derivatives for a dihedral for a given angle. */
  Utils::AutomaticDifferentiation::Second1D getInteraction(double angle) const;
  /**
   * @brief Returns whether parameters are available for this bond.
   */
  bool hasParameters() const;
  /** @brief Setter for the cosine pre-factor (typically: -1 for SFAM, +1 for GAFF) */
  void setCosinePreFactor(double cosPreFactor);

 private:
  double halfBarrierHeight_;
  int periodicity_;
  double phaseShift_;
  bool parametersAreAvailable_;
  double cosPreFactor_ = -1.0; // Correct for SFAM, but has to be set to +1.0 for GAFF.
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_DIHEDRAL_H