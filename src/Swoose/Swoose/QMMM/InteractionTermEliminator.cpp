/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "InteractionTermEliminator.h"
#include <Swoose/MolecularMechanics/GAFF/GaffMolecularMechanicsCalculator.h>
#include <Swoose/MolecularMechanics/SFAM/SfamMolecularMechanicsCalculator.h>

namespace Scine {
namespace Qmmm {

InteractionTermEliminator::InteractionTermEliminator(const std::vector<int>& listOfQmAtoms,
                                                     std::shared_ptr<MolecularMechanics::MolecularMechanicsCalculator> calculator)
  : qmAtoms_(listOfQmAtoms.begin(), listOfQmAtoms.end()), calculator_(calculator) {
}

void InteractionTermEliminator::eliminateInteractionTerms(bool electrostaticEmbedding, bool eliminateEnvironmentOnlyTerms) {
  eliminateEnvironmentOnlyTerms_ = eliminateEnvironmentOnlyTerms;
  auto calculatorName = calculator_->name();
  if (calculatorName == "SFAM") {
    auto calc = std::dynamic_pointer_cast<MolecularMechanics::SfamMolecularMechanicsCalculator>(calculator_);
    if (calc == nullptr)
      throw std::runtime_error("SFAM Calculator could not be casted to derived class.");
    eliminateSharedInteractionTerms(electrostaticEmbedding, *calc);
    eliminateImproperDihedralTerms(calc->improperDihedralsEvaluator_->improperDihedrals_);
    eliminateDispersionTerms(calc->dispersionEvaluator_->dispersions_);
    eliminateRepulsionTerms(calc->repulsionEvaluator_->repulsions_);
    eliminateHydrogenBondTerms(calc->hydrogenBondEvaluator_->hydrogenBondTerms_);
  }
  else if (calculatorName == "GAFF") {
    auto calc = std::dynamic_pointer_cast<MolecularMechanics::GaffMolecularMechanicsCalculator>(calculator_);
    if (calc == nullptr)
      throw std::runtime_error("GAFF Calculator could not be casted to derived class.");
    eliminateSharedInteractionTerms(electrostaticEmbedding, *calc);
    eliminateDihedralTerms(calc->improperDihedralsEvaluator_->dihedrals_, true);
    eliminateLennardJonesTerms(calc->lennardJonesEvaluator_->ljTerms_);
  }
  else {
    throw std::runtime_error("The given MM model is not supported by the interaction term eliminator.");
  }
}

template<typename CalculatorType>
void InteractionTermEliminator::eliminateSharedInteractionTerms(bool electrostaticEmbedding, CalculatorType& calculator) {
  eliminateBondedTerms(calculator.bondsEvaluator_->bonds_);
  eliminateAngleTerms(calculator.anglesEvaluator_->angles_);
  eliminateDihedralTerms(calculator.dihedralsEvaluator_->dihedrals_);
  eliminateElectrostaticTerms(calculator.electrostaticEvaluator_->electrostaticTerms_, electrostaticEmbedding);
}

void InteractionTermEliminator::eliminateTerm(MolecularMechanics::InteractionTermBase& term,
                                              const std::vector<int>& atomsInTerm, int allowedInQmRegion) {
  int inQmRegion = 0;
  for (const auto& atom : atomsInTerm) {
    if (isQmAtom(atom))
      inQmRegion++;
  }
  if (inQmRegion > allowedInQmRegion)
    term.disable();
  if (eliminateEnvironmentOnlyTerms_) {
    if (inQmRegion == 0)
      term.disable();
  }
}

bool InteractionTermEliminator::isQmAtom(int index) {
  return qmAtoms_.find(index) != qmAtoms_.end();
}

void InteractionTermEliminator::eliminateBondedTerms(std::vector<MolecularMechanics::BondedTerm>& bondedTerms) {
  for (auto& bondedTerm : bondedTerms) {
    std::vector<int> atoms = {bondedTerm.getFirstAtom(), bondedTerm.getSecondAtom()};
    if (eliminateEnvironmentOnlyTerms_)
      bondedTerm.disable();
    else
      eliminateTerm(bondedTerm, atoms, 1);
  }
}

void InteractionTermEliminator::eliminateAngleTerms(std::vector<MolecularMechanics::AngleTerm>& angleTerms) {
  for (auto& angleTerm : angleTerms) {
    std::vector<int> atoms = {angleTerm.getFirstAtom(), angleTerm.getSecondAtom(), angleTerm.getThirdAtom()};
    if (eliminateEnvironmentOnlyTerms_)
      angleTerm.disable();
    else
      eliminateTerm(angleTerm, atoms, 1); // TODO: 1 or 2?
  }
}

void InteractionTermEliminator::eliminateDihedralTerms(std::vector<MolecularMechanics::DihedralTerm>& dihedralTerms,
                                                       bool gaffImproper) {
  for (auto& dihedralTerm : dihedralTerms) {
    std::vector<int> atoms = {dihedralTerm.getFirstAtom(), dihedralTerm.getSecondAtom(), dihedralTerm.getThirdAtom(),
                              dihedralTerm.getFourthAtom()};
    if (eliminateEnvironmentOnlyTerms_)
      dihedralTerm.disable();
    else {
      int allowedInQmRegion = gaffImproper ? 0 : 2; // TODO: 2 or 3?
      eliminateTerm(dihedralTerm, atoms, allowedInQmRegion);
    }
  }
}

void InteractionTermEliminator::eliminateImproperDihedralTerms(std::vector<MolecularMechanics::ImproperDihedralTerm>& improperDihedralTerms) {
  for (auto& improperDihedralTerm : improperDihedralTerms) {
    std::vector<int> atoms = {improperDihedralTerm.getCentralAtom()};
    if (eliminateEnvironmentOnlyTerms_)
      improperDihedralTerm.disable();
    else
      eliminateTerm(improperDihedralTerm, atoms, 0);
  }
}

void InteractionTermEliminator::eliminateDispersionTerms(std::vector<MolecularMechanics::DispersionTerm>& dispersionTerms) {
  for (auto& dispersionTerm : dispersionTerms) {
    std::vector<int> atoms = {dispersionTerm.getFirstAtom(), dispersionTerm.getSecondAtom()};
    eliminateTerm(dispersionTerm, atoms, 1);
  }
}

void InteractionTermEliminator::eliminateRepulsionTerms(std::vector<MolecularMechanics::RepulsionTerm>& repulsionTerms) {
  for (auto& repulsionTerm : repulsionTerms) {
    std::vector<int> atoms = {repulsionTerm.getFirstAtom(), repulsionTerm.getSecondAtom()};
    eliminateTerm(repulsionTerm, atoms, 1);
  }
}

void InteractionTermEliminator::eliminateLennardJonesTerms(std::vector<MolecularMechanics::LennardJonesTerm>& ljTerms) {
  for (auto& ljTerm : ljTerms) {
    std::vector<int> atoms = {ljTerm.getFirstAtom(), ljTerm.getSecondAtom()};
    eliminateTerm(ljTerm, atoms, 1);
  }
}

void InteractionTermEliminator::eliminateHydrogenBondTerms(std::vector<MolecularMechanics::HydrogenBondTerm>& hydrogenBondTerms) {
  for (auto& hydrogenBondTerm : hydrogenBondTerms) {
    std::vector<int> atoms = {hydrogenBondTerm.getHydrogenAtom(), hydrogenBondTerm.getDonorAtom(),
                              hydrogenBondTerm.getAcceptorAtom()};
    if (eliminateEnvironmentOnlyTerms_)
      hydrogenBondTerm.disable();
    else
      eliminateTerm(hydrogenBondTerm, atoms, 2); // Eliminate only if all three atoms are inside the QM region
  }
}

void InteractionTermEliminator::eliminateElectrostaticTerms(std::vector<MolecularMechanics::ElectrostaticTerm>& electrostaticTerms,
                                                            bool electrostaticEmbedding) {
  for (auto& electrostaticTerm : electrostaticTerms) {
    std::vector<int> atoms = {electrostaticTerm.getFirstAtom(), electrostaticTerm.getSecondAtom()};
    // If electrostatic embedding is switched on:
    // Eliminate if at least one atom is in the QM region,
    // because QM-MM electrostatic interaction is covered by the electrostatic embedding
    if (electrostaticEmbedding)
      eliminateTerm(electrostaticTerm, atoms, 0);
    else // For mechanical embedding: If one atom is in the QM region -> don't eliminate this term
      eliminateTerm(electrostaticTerm, atoms, 1);
  }
}

void InteractionTermEliminator::reset() {
  auto calculatorName = calculator_->name();
  if (calculatorName == "SFAM") {
    auto calc = std::dynamic_pointer_cast<MolecularMechanics::SfamMolecularMechanicsCalculator>(calculator_);
    if (calc == nullptr)
      throw std::runtime_error("SFAM Calculator could not be casted to derived class.");
    enableSharedInteractionTerms(*calc);
    enableTerms(calc->improperDihedralsEvaluator_->improperDihedrals_);
    enableTerms(calc->dispersionEvaluator_->dispersions_);
    enableTerms(calc->repulsionEvaluator_->repulsions_);
    enableTerms(calc->hydrogenBondEvaluator_->hydrogenBondTerms_);
  }
  else if (calculatorName == "GAFF") {
    auto calc = std::dynamic_pointer_cast<MolecularMechanics::GaffMolecularMechanicsCalculator>(calculator_);
    if (calc == nullptr)
      throw std::runtime_error("GAFF Calculator could not be casted to derived class.");
    enableSharedInteractionTerms(*calc);
    enableTerms(calc->improperDihedralsEvaluator_->dihedrals_);
    enableTerms(calc->lennardJonesEvaluator_->ljTerms_);
  }
  else {
    throw std::runtime_error("The given MM model is not supported by the interaction term eliminator.");
  }
}

template<typename T>
void InteractionTermEliminator::enableTerms(std::vector<T>& terms) {
  for (auto& t : terms) {
    t.enable();
  }
}

template<typename CalculatorType>
void InteractionTermEliminator::enableSharedInteractionTerms(CalculatorType& calculator) {
  enableTerms(calculator.bondsEvaluator_->bonds_);
  enableTerms(calculator.anglesEvaluator_->angles_);
  enableTerms(calculator.dihedralsEvaluator_->dihedrals_);
  enableTerms(calculator.electrostaticEvaluator_->electrostaticTerms_);
}

} // namespace Qmmm
} // namespace Scine
