/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_DISPERSIONC6PARAMETERS_H
#define MOLECULARMECHANICS_DISPERSIONC6PARAMETERS_H

namespace Scine {

namespace Utils {
class AtomCollection;
} // namespace Utils

namespace MolecularMechanics {
class SfamParameters;
class AtomTypesHolder;

namespace DispersionC6Parameters {

/**
 * @brief Fills the C6 parameter matrix in the MM parameters with values calculated for the current structure.
 * @param parameters The MM parameters to manipulate.
 * @param structure The current molecular structure.
 * @param atomTypes The atom types corresponding to the molecular structure.
 */
void fillC6MatrixForCurrentStructure(SfamParameters& parameters, const Utils::AtomCollection& structure,
                                     const AtomTypesHolder& atomTypes);

} // namespace DispersionC6Parameters
} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_DISPERSIONC6PARAMETERS_H
