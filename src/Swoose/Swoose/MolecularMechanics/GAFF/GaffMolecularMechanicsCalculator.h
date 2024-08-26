/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef GAFFMOLECULARMECHANICSCALCULATOR_H
#define GAFFMOLECULARMECHANICSCALCULATOR_H

#include "../Interactions/AnglesEvaluator.h"
#include "../Interactions/BondsEvaluator.h"
#include "../Interactions/DihedralsEvaluator.h"
#include "../Interactions/ElectrostaticEvaluator.h"
#include "../Interactions/LennardJonesEvaluator.h"
#include "../MolecularMechanicsCalculator.h"
#include "../Topology/IndexedStructuralTopology.h"
#include "GaffAtomTypeIdentifier.h"
#include "GaffParameters.h"

namespace Scine {

namespace Qmmm {
class InteractionTermEliminator;
} // namespace Qmmm

namespace MolecularMechanics {
class GaffParameters;
class AtomTypesHolder;
struct BondType;
struct AngleType;
struct DihedralType;
struct ImproperDihedralType;
/**
 * @class GaffMolecularMechanicsCalculator GaffMolecularMechanicsCalculator.h
 * @brief Calculator for the GAFF Molecular Mechanics method.
 *
 *        Note: First call the constructor, then set parameter and connectivity files in the settings
 *        and after that one can call "setStructure", because "setStructure" does the initialization,
 *        which reads in the parameters and connectivity from the previously specified files.
 */
class GaffMolecularMechanicsCalculator final
  : public Utils::CloneInterface<GaffMolecularMechanicsCalculator, MolecularMechanicsCalculator, Core::Calculator> {
 public:
  static constexpr const char* model = "GAFF";
  /// @brief Constructor.
  GaffMolecularMechanicsCalculator();
  /// @brief Default Destructor.
  ~GaffMolecularMechanicsCalculator() override = default;
  /// @brief Copy Constructor.
  GaffMolecularMechanicsCalculator(const GaffMolecularMechanicsCalculator& rhs);
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
  friend class Qmmm::InteractionTermEliminator;
  // Setter for the internal copy of the parameters.
  void setParameters(GaffParameters parameters);
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
  std::unique_ptr<DihedralsEvaluator> improperDihedralsEvaluator_;
  std::unique_ptr<ElectrostaticEvaluator> electrostaticEvaluator_;
  std::unique_ptr<LennardJonesEvaluator> lennardJonesEvaluator_;
  bool printContributionsMolecularMechanics_;
  std::string atomicChargesFile_;
  std::string atomTypesFile_;
  bool onlyCalculateBondedContribution_;
  bool detectBondsWithCovalentRadii_;
  /*
   * Decides whether the potential terms that correspond to non-bonded interaction beyond the cutoff radius,
   * are ignored during the initialization stage of the calculator. This can result in a great speed-up, however,
   * the calculator may has to be re-initialized when the structure changes significantly.
   */
  bool applyCutoffDuringInitialization_;
  double nonCovalentCutoffRadius_;
  std::string connectivityFilePath_;
  std::string parameterFilePath_;
  /*
   * Tracks whether the parameter file path has been recently changed.
   * This boolean is set to true when the member 'parameterFilePath_' is updated
   * and it is set back to false when these parameters are parsed and stored in the
   * 'parameters_' member.
   */
  bool parameterFilePathHasBeenChanged_ = false;
  // Whether this calculator is currently calculating a Hessian -> no detailed output printing
  bool hessianMode_ = false;
  // Copy of the parameters as a member that is also cloned when the MM calculator is cloned.
  GaffParameters parameters_;
  bool parametersHaveBeenSetInternally_ = false;
  // Private methods:
  void initialize();
  void generatePotentialTerms(const GaffParameters& parameters, const IndexedStructuralTopology& topology,
                              const AtomTypesHolder& atomTypes);
  void generatePotentialTerms(const std::string& parameterPath);
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // GAFFMOLECULARMECHANICSCALCULATOR_H
