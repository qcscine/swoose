/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_POTENTIALTERMSHELPER_H
#define MOLECULARMECHANICS_POTENTIALTERMSHELPER_H

#include "Interactions/AngleTerm.h"
#include "Interactions/BondedTerm.h"
#include <Eigen/Dense>
#include <vector>

namespace Scine {
namespace MolecularMechanics {
class IndexedStructuralTopology;
class AtomTypesHolder;
class MMParameters;

namespace PotentialTermsHelper {

/**
 * @brief Calculates and returns the exclusion-type matrix, which encodes for each atom pair whether its
 *        non-covalent contribution is included, scaled or excluded. The matrix holds integers that encode
 *        the following meaning: 0 = excluded, 1 = included, -1 = scaled.
 * @param topology The topology of the system.
 * @param nAtoms The number of atoms of the system.
 * @return Exclusion-type matrix, that encodes for each atom pair how its non-covalent contribution is handled.
 */
Eigen::MatrixXi getExclusionTypeMatrix(const IndexedStructuralTopology& topology, int nAtoms);

/**
 * @brief Getter for the vector of bonded terms. This is the same for GAFF and SFAM.
 * @param topology The topology of the system.
 * @param parameters The MM parameters container.
 * @param atomTypesHolder The atom types of the system.
 * @return Vector of bonded terms.
 */
std::vector<BondedTerm> getBondedTerms(const IndexedStructuralTopology& topology, const MMParameters& parameters,
                                       const AtomTypesHolder& atomTypesHolder);

/**
 * @brief Getter for the vector of angle terms. This is the same for GAFF and SFAM.
 * @param topology The topology of the system.
 * @param parameters The MM parameters container.
 * @param atomTypesHolder The atom types of the system.
 * @return Vector of angle terms.
 */
std::vector<AngleTerm> getAngleTerms(const IndexedStructuralTopology& topology, const MMParameters& parameters,
                                     const AtomTypesHolder& atomTypesHolder);
} // namespace PotentialTermsHelper
} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_POTENTIALTERMSHELPER_H
