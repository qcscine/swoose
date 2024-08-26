/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_GAFFPOTENTIALTERMSGENERATOR_H
#define MOLECULARMECHANICS_GAFFPOTENTIALTERMSGENERATOR_H

#include "../Interactions/AngleTerm.h"
#include "../Interactions/BondedTerm.h"
#include "../Interactions/DihedralTerm.h"
#include "../Interactions/ImproperDihedralTerm.h"
#include <vector>

namespace Scine {

namespace Core {
class Log;
} // namespace Core

namespace MolecularMechanics {
class IndexedStructuralTopology;
class AtomTypesHolder;
class GaffParameters;

/**
 * @class GaffPotentialTermsGenerator GaffPotentialTermsGenerator.h
 * @brief This class creates the actual potential terms that are employed
 *        during a GAFF molecular mechanics calculation.
 *        It combines information about structure (bonds, angles, etc),
 *        atom types and parameters in order to do so.
 */
class GaffPotentialTermsGenerator {
 public:
  /**
   * @brief Constructor.
   * @param nAtoms Number of atoms in the system.
   * @param atomTypes The atom types.
   * @param topology The topology.
   * @param parameters The parameters of GAFF.
   * @param log The logger.
   */
  GaffPotentialTermsGenerator(int nAtoms, const AtomTypesHolder& atomTypes, const IndexedStructuralTopology& topology,
                              const GaffParameters& parameters, const Utils::PositionCollection& positions, Core::Log& log);

  /** @brief Getter for bonded terms */
  std::vector<BondedTerm> getBondedTerms();
  /** @brief Getter for angle terms */
  std::vector<AngleTerm> getAngleTerms();
  /** @brief Getter for dihedral terms */
  std::vector<DihedralTerm> getDihedralTerms();
  /** @brief Getter for improper dihedral terms. In GAFF, an improper dihedral is treated as a normal dihedral. */
  std::vector<DihedralTerm> getImproperDihedralTerms();

 private:
  int nAtoms_;
  const AtomTypesHolder& atomTypesHolder_;
  const IndexedStructuralTopology& topology_;
  const GaffParameters& parameters_;
  const Utils::PositionCollection& positions_;
  // The logger
  Core::Log& log_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_GAFFPOTENTIALTERMSGENERATOR_H
