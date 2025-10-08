/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_LENNARDJONESPARAMETERS_H
#define MOLECULARMECHANICS_LENNARDJONESPARAMETERS_H

#include "../Interactions/LennardJones.h"

namespace Scine {
namespace MolecularMechanics {

/**
 * @class LennardJonesParameters LennardJonesParameters.h
 * @brief Class containing the parameters for a van der Waals LJ-type interaction.
 */
class LennardJonesParameters {
 public:
  /**
   * @brief Constructor.
   * @param vdwRadius Unit: Angstrom
   * @param wellDepth Unit: kcal/mol
   */
  LennardJonesParameters(double vdwRadius, double wellDepth);
  /**
   * @brief Getter for the van-der-Waals radius in bohr.
   * @return The van-der-Waals radius in bohr.
   */
  double getVdwRadius() const;
  /**
   * @brief Getter for the well depth in hartree.
   * @return The potential well depth.
   */
  double getWellDepth() const;
  /**
   * @brief Getter for the atom type index.
   *        The index in the atom type list can be cached by this class.
   * @return The atom type index.
   */
  unsigned int getAtomTypeIndex() const;
  /**
   * @brief Setter for the atom type index.
   *        The index in the atom type list can be cached by this class.
   * @param index The index.
   */
  void setAtomTypeIndex(unsigned int index);

 private:
  double vdwRadius_; // Unit: Bohr
  double wellDepth_; // Unit: Hartree
  unsigned int atomTypeIndex_ = 0;
};

class LennardJonesPairParameters {
 public:
  LennardJonesPairParameters(const LennardJonesParameters& first, const LennardJonesParameters& second);

  const double& getACoeff();
  const double& getBCoeff();
  Utils::AutomaticDifferentiation::Second1D getInteraction(const double& distance, const double& scaling) const;

 private:
  double aCoeff_;
  double bCoeff_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_BONDPARAMETERS_H
