/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_INTERACTIONTERMELIMINATOR_H
#define SWOOSE_QMMM_INTERACTIONTERMELIMINATOR_H

#include <memory>
#include <unordered_set>
#include <vector>

namespace Scine {

namespace MolecularMechanics {
class MolecularMechanicsCalculator;
class BondedTerm;
class AngleTerm;
class DihedralTerm;
class ImproperDihedralTerm;
class DispersionTerm;
class RepulsionTerm;
class ElectrostaticTerm;
class HydrogenBondTerm;
class LennardJonesTerm;
class InteractionTermBase;
} // namespace MolecularMechanics

namespace Qmmm {

/**
 * @class InteractionTermEliminator InteractionTermEliminator.h
 * @brief This class handles the elimination of MM interaction terms, which are already covered
 *        by the QM calculation in QM/MM.
 */
class InteractionTermEliminator {
 public:
  /// @brief Constructor.
  InteractionTermEliminator(const std::vector<int>& listOfQmAtoms,
                            std::shared_ptr<MolecularMechanics::MolecularMechanicsCalculator> calculator);

  /**
   * @brief Eliminates the interactions terms in an MM calculator for QM/MM,
   *        which are already covered by the QM calculation.
   * @param electrostaticEmbedding Whether electrostatic embedding is switched on in the QM/MM calculator.
   * @param eliminateEnvironmentOnlyTerms Whether the terms that are only
   *                                      located within the environment shall be eliminated.
   */
  void eliminateInteractionTerms(bool electrostaticEmbedding, bool eliminateEnvironmentOnlyTerms = false);

  /**
   * @brief Enables all of the interactions terms (reverting elimination and therefore resetting the state of the MM).
   */
  void reset();

 private:
  // @brief Eliminates those bonded interactions which are already covered by the QM calculation in QM/MM.
  void eliminateBondedTerms(std::vector<MolecularMechanics::BondedTerm>& bondedTerms);

  // @brief Eliminates those angle interactions which are already covered by the QM calculation in QM/MM.
  void eliminateAngleTerms(std::vector<MolecularMechanics::AngleTerm>& angleTerms);

  /*
   * @brief Eliminates those dihedral interactions which are already covered by the QM calculation in QM/MM.
   *        The additional boolean argument "gaffImproper" is needed if this function is called for the GAFF
   *        method to eliminate improper dihedrals (which are treated as normal dihedrals in GAFF). In such a case,
   *        the allowed number of QM atoms is different.
   */
  void eliminateDihedralTerms(std::vector<MolecularMechanics::DihedralTerm>& dihedralTerms, bool gaffImproper = false);

  // @brief Eliminates those improper dihedral interactions which are already covered by the QM calculation in QM/MM.
  void eliminateImproperDihedralTerms(std::vector<MolecularMechanics::ImproperDihedralTerm>& improperDihedralTerms);

  // @brief Eliminates those dispersion interactions which are already covered by the QM calculation in QM/MM.
  void eliminateDispersionTerms(std::vector<MolecularMechanics::DispersionTerm>& dispersionTerms);

  // @brief Eliminates those Pauli repulsion interactions which are already covered by the QM calculation in QM/MM.
  void eliminateRepulsionTerms(std::vector<MolecularMechanics::RepulsionTerm>& repulsionTerms);

  // @brief Eliminates those Lennard-Jones interactions which are already covered by the QM calculation in QM/MM.
  void eliminateLennardJonesTerms(std::vector<MolecularMechanics::LennardJonesTerm>& ljTerms);

  // @brief Eliminates those hydrogen bond interactions which are already covered by the QM calculation in QM/MM.
  void eliminateHydrogenBondTerms(std::vector<MolecularMechanics::HydrogenBondTerm>& hydrogenBondTerms);

  /*
   * @brief Eliminates those electrostatic interactions which are already covered by the QM calculation in QM/MM.
   * @param electrostaticTerms The electrostatic terms of the MM model.
   * @param electrostaticEmbedding Whether electrostatic embedding is switched on in the QM/MM calculator.
   */
  void eliminateElectrostaticTerms(std::vector<MolecularMechanics::ElectrostaticTerm>& electrostaticTerms,
                                   bool electrostaticEmbedding);

  /*
   * @brief Eliminates the given interaction term if more than 'allowedInQmRegion' of its atoms are in the QM region.
   */
  void eliminateTerm(MolecularMechanics::InteractionTermBase& term, const std::vector<int>& atomsInTerm, int allowedInQmRegion);

  /*
   * @brief Returns whether an atom with the given index is in the QM region.
   */
  bool isQmAtom(int index);

  // Enables all representatives of one type of interaction terms.
  template<typename T>
  void enableTerms(std::vector<T>& terms);

  // Eliminates all of the shared interaction terms that are present in all of the possible MM calculators.
  template<typename CalculatorType>
  void eliminateSharedInteractionTerms(bool electrostaticEmbedding, CalculatorType& calculator);

  // Enables all of the shared interaction terms that are present in all of the possible MM calculators.
  template<typename CalculatorType>
  void enableSharedInteractionTerms(CalculatorType& calculator);

  // Hash table containing the indices of the atoms in the QM region.
  std::unordered_set<int> qmAtoms_;

  // The MM calculator
  std::shared_ptr<MolecularMechanics::MolecularMechanicsCalculator> calculator_;

  // Boolean to keep track of whether we are in a reduced QM/MM energy calculation.
  bool eliminateEnvironmentOnlyTerms_ = false;
};

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_INTERACTIONTERMELIMINATOR_H
