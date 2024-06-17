/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_SFAMPOTENTIALTERMSGENERATOR_H
#define MOLECULARMECHANICS_SFAMPOTENTIALTERMSGENERATOR_H

#include "../Interactions/AngleTerm.h"
#include "../Interactions/BondedTerm.h"
#include "../Interactions/DihedralTerm.h"
#include "../Interactions/DispersionTerm.h"
#include "../Interactions/ElectrostaticTerm.h"
#include "../Interactions/HydrogenBondTerm.h"
#include "../Interactions/ImproperDihedralTerm.h"
#include "../Interactions/RepulsionTerm.h"
#include <vector>

namespace Scine {

namespace Core {
class Log;
} // namespace Core

namespace MolecularMechanics {
class IndexedStructuralTopology;
class AtomTypesHolder;
class SfamParameters;

/**
 * @class SfamPotentialTermsGenerator SfamPotentialTermsGenerator.h
 * @brief This class creates the actual potential terms that are employed
 *        during an SFAM molecular mechanics calculation.
 *        It combines information about structure (bonds, angles, etc),
 *        atom types and parameters in order to do so.
 */
class SfamPotentialTermsGenerator {
 public:
  /**
   * @brief Constructor.
   * @param nAtoms Number of atoms in the system.
   * @param atomTypes The atom types.
   * @param topology The topology.
   * @param parameters The parameters of SFAM's MM model.
   * @param nonCovalentCutoffRadius The cutoff radius for non-covalent interactions.
   * @param log The logger.
   */
  SfamPotentialTermsGenerator(int nAtoms, const AtomTypesHolder& atomTypes, const IndexedStructuralTopology& topology,
                              const SfamParameters& parameters, const Utils::PositionCollection& positions,
                              const double& nonCovalentCutoffRadius, Core::Log& log);

  /** @brief Getter for bonded terms */
  std::vector<BondedTerm> getBondedTerms();
  /** @brief Getter for angle terms */
  std::vector<AngleTerm> getAngleTerms();
  /** @brief Getter for dihedral terms */
  std::vector<DihedralTerm> getDihedralTerms();
  /** @brief Getter for improper dihedral terms */
  std::vector<ImproperDihedralTerm> getImproperDihedralTerms();
  /** @brief Getter for dispersion terms with decision whether to use a cutoff radius. */
  std::vector<DispersionTerm> getDispersionTerms(bool applyCutoff);
  /** @brief Getter for Pauli repulsion terms with decision whether to use a cutoff radius. */
  std::vector<RepulsionTerm> getRepulsionTerms(bool applyCutoff);
  /** @brief Getter for electrostatic terms with decision whether to use a cutoff radius. */
  std::vector<ElectrostaticTerm> getElectrostaticTerms(bool applyCutoff);
  /** @brief Getter for hydrogen bond terms */
  std::vector<HydrogenBondTerm> getHydrogenBondTerms();

 private:
  int nAtoms_;
  const AtomTypesHolder& atomTypesHolder_;
  const IndexedStructuralTopology& topology_;
  const SfamParameters& parameters_;
  const Utils::PositionCollection& positions_;
  /*
   * @brief The cutoff radius for the non-covalent interactions.
   *
   *        It is a pointer so that it is updated with the settings.
   *        This is important since it is passed on from here to the actual DispersionTerms, etc.
   *        A reference does not work because then the assignment operator of DispersionTerms, etc. vanishes.
   */
  std::shared_ptr<double> cutoff_;
  // The logger
  Core::Log& log_;

  static constexpr double scalingFactorForElectrostaticOneFourTerms_ = 1.0;
  static constexpr double scalingFactorForDispersionOneFourTerms_ = 1.0;
  static constexpr double scalingFactorForRepulsionOneFourTerms_ = 1.0;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_SFAMPOTENTIALTERMSGENERATOR_H
