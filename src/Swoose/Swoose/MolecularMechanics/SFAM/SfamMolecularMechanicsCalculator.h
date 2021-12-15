/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SFAMMOLECULARMECHANICSCALCULATOR_H
#define SFAMMOLECULARMECHANICSCALCULATOR_H

#include "../Interactions/AnglesEvaluator.h"
#include "../Interactions/BondsEvaluator.h"
#include "../Interactions/DihedralsEvaluator.h"
#include "../Interactions/DispersionEvaluator.h"
#include "../Interactions/ElectrostaticEvaluator.h"
#include "../Interactions/HydrogenBondEvaluator.h"
#include "../Interactions/ImproperDihedralsEvaluator.h"
#include "../Interactions/RepulsionEvaluator.h"
#include "../MolecularMechanicsCalculator.h"
#include "../Topology/IndexedStructuralTopology.h"
#include "SfamAtomTypeIdentifier.h"
#include "SfamParameters.h"

namespace Scine {

namespace MMParametrization {
class UpdateFunctionManager;
} // namespace MMParametrization

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace MolecularMechanics {
class SfamParameters;
class AtomTypesHolder;
class BondType;
class AngleType;
class DihedralType;
class ImproperDihedralType;
/**
 * @class SfamMolecularMechanicsCalculator SfamMolecularMechanicsCalculator.h
 * @brief Calculator for the SFAM Molecular Mechanics method.
 *
 *        Note: First call the constructor, then set parameter and connectivity files in the settings
 *        and after that one can call "setStructure", because "setStructure" does the initialization,
 *        which reads in the parameters and connectivity from the previously specified files.
 */
class SfamMolecularMechanicsCalculator final
  : public Utils::CloneInterface<SfamMolecularMechanicsCalculator, MolecularMechanicsCalculator> {
 public:
  static constexpr const char* model = "SFAM";
  /// @brief Constructor.
  SfamMolecularMechanicsCalculator();
  /// @brief Default Destructor.
  ~SfamMolecularMechanicsCalculator() override = default;
  /// @brief Copy Constructor.
  SfamMolecularMechanicsCalculator(const SfamMolecularMechanicsCalculator& rhs);
  /**
   * @brief Changes the molecular structure to calculate.
   * @param structure A new Utils::AtomCollection to save.
   */
  void setStructure(const Utils::AtomCollection& structure) override;
  /**
   * @brief The main function running calculations.
   *
   * @param description   The calculation description.
   * @return Utils::Result Return the result of the calculation. The object contains the
   *                       properties that were given as requirement by the
   *                       Calculator::setRequiredProperties function.
   */
  const Utils::Results& calculate(std::string description) override;
  /**
   * @brief Getter for the name of the Calculator.
   * @return Returns the name of the Calculator.
   */
  std::string name() const override;

 private:
  // friend class declarations
  friend class MMParametrization::UpdateFunctionManager;
  friend class Qmmm::InteractionTermEliminator;
  // This function is only used by the friend class in the MM parametrization algorithm.
  void setListsOfNeighbors(std::vector<std::list<int>> listsOfNeighbors);
  // Setter for the internal copy of the parameters.
  void setParameters(SfamParameters parameters);
  /*
   * @brief Implementation of a calculation.
   */
  const Utils::Results& calculateImpl(std::string description);
  /*
   * @brief Apply settings.
   */
  void applySettings();

  // Other private members:
  std::unique_ptr<BondsEvaluator> bondsEvaluator_;
  std::unique_ptr<AnglesEvaluator> anglesEvaluator_;
  std::unique_ptr<DihedralsEvaluator> dihedralsEvaluator_;
  std::unique_ptr<ImproperDihedralsEvaluator> improperDihedralsEvaluator_;
  std::unique_ptr<DispersionEvaluator> dispersionEvaluator_;
  std::unique_ptr<RepulsionEvaluator> repulsionEvaluator_;
  std::unique_ptr<ElectrostaticEvaluator> electrostaticEvaluator_;
  std::unique_ptr<HydrogenBondEvaluator> hydrogenBondEvaluator_;
  bool printContributionsMolecularMechanics_;
  bool onlyCalculateBondedContribution_;
  bool detectBondsWithCovalentRadii_;
  bool includeHydrogenBonds_;
  /*
   * Decides whether the potential terms that correspond to non-bonded interaction beyond the cutoff radius,
   * are ignored during the initialization stage of the calculator. This can result in a great speed-up, however,
   * the calculator may has to be re-initialized when the structure changes significantly.
   */
  bool applyCutoffDuringInitialization_;
  double nonCovalentCutoffRadius_;
  SfamAtomTypeLevel sfamAtomTypeLevel_;
  std::string connectivityFilePath_;
  std::string parameterFilePath_;
  /*
   * Tracks whether the parameter file path has been recently changed.
   * This boolean is set to true when the member 'parameterFilePath_' is updated
   * and it is set back to false when these parameters are parsed and stored in the
   * 'parameters_' member.
   */
  bool parameterFilePathHasBeenChanged_ = false;
  // Copy of the parameters as a member that is also cloned when the MM calculator is cloned.
  SfamParameters parameters_;
  bool parametersHaveBeenSetInternally_ = false;
  // Whether this calculator is currently calculating a Hessian -> no detailed output printing
  bool hessianMode_ = false;
  // Private methods:
  void initialize();
  void generatePotentialTerms(const SfamParameters& parameters, const IndexedStructuralTopology& topology,
                              const AtomTypesHolder& atomTypes);
  void generatePotentialTerms(const std::string& parameterPath);
  /*
   * A list of atom indices for which the contributions in the numerical part of the Hessian
   * calculations are considered. All other atom contributions in this Hessian will be ignored.
   * This is useful for speeding up the parametrizations of the model.
   */
  std::vector<int> atomsToConsiderForHessian_ = {};
  // Setter for the 'atomsToConsiderForHessian_' member
  void setAtomsToConsiderForHessian(std::vector<int> atomsToConsiderForHessian);
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // SFAMMOLECULARMECHANICSCALCULATOR_H
