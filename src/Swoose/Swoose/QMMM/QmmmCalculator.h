/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMMMCALCULATOR_H
#define SWOOSE_QMMM_QMMMCALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/EmbeddingCalculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Technical/CloneInterface.h>
#include <algorithm>

namespace Scine {

namespace MolecularMechanics {
class MolecularMechanicsCalculator;
} // namespace MolecularMechanics

namespace Qmmm {

/**
 * @class QmmmCalculator QmmmCalculator.h
 * @brief Calculator implementing the QM/MM method.
 *
 */
class QmmmCalculator : public Utils::CloneInterface<QmmmCalculator, Core::EmbeddingCalculator, Core::Calculator> {
 public:
  static constexpr const char* model = "QM-SFAM";
  /// @brief Constructor.
  QmmmCalculator();
  /// @brief Destructor.
  ~QmmmCalculator() override;
  /// @brief Copy Constructor.
  QmmmCalculator(const QmmmCalculator& rhs);
  /**
   * @brief Sets the underlying QM and MM calculators.
   */
  void setUnderlyingCalculators(std::vector<std::shared_ptr<Core::Calculator>> underlyingCalculators) override;
  /**
   * @brief Sets the underlying QM and MM calculators.
   */
  void addUnderlyingSettings() override;
  /**
   * @brief Gets the underlying QM and MM calculators.
   */
  std::vector<std::shared_ptr<Core::Calculator>> getUnderlyingCalculators() const override;
  /**
   * @brief Changes the molecular structure to calculate.
   * @param structure A new Utils::AtomCollection to save.
   */
  void setStructure(const Utils::AtomCollection& structure) override;
  /**
   * @brief Gets the molecular structure as a const Utils::AtomCollection&.
   * @return a const Utils::AtomCollection&.
   */
  std::unique_ptr<Utils::AtomCollection> getStructure() const override;
  /**
   * @brief Allows to modify the positions of the underlying Utils::AtomCollection.
   * @param newPositions The new positions to be assigned to the underlying Utils::AtomCollection.
   */
  void modifyPositions(Utils::PositionCollection newPositions) override;
  /**
   * @brief Getter for the coordinates of the underlying Utils::AtomCollection.
   */
  const Utils::PositionCollection& getPositions() const override;
  /**
   * @brief Sets the properties to calculate.
   * @param requiredProperties A Utils::PropertyList, a sequence of bits that represent the
   *        properties that must be calculated.
   */
  void setRequiredProperties(const Utils::PropertyList& requiredProperties) override;
  /**
   * @brief Getter for the properties to calculate.
   */
  Utils::PropertyList getRequiredProperties() const override;
  /**
   * @brief Returns the list of the possible properties to calculate.
   */
  Utils::PropertyList possibleProperties() const override;
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
  /**
   * @brief Accessor for the settings.
   * @return Utils::Settings& The settings.
   */
  Utils::Settings& settings() override;
  /**
   * @brief Constant accessor for the settings.
   * @return const Utils::Settings& The settings.
   */
  const Utils::Settings& settings() const override;
  /**
   * @brief Accessor for the saved instance of Utils::Results.
   * @return Utils::Results& The results of the previous calculation.
   */
  Utils::Results& results() override;
  /**
   * @brief Constant accessor for the Utils::Results.
   * @return const Utils::Results& The results of the previous calculation.
   */
  const Utils::Results& results() const override;
  /**
   * @brief Whether the calculator supports a method family.
   * @param methodFamily Identifier for the method family.
   * @return Whether the calculator supports a method family.
   */
  bool supportsMethodFamily(const std::string& methodFamily) const override;
  /**
   * @brief Implements Core::StateHandableObject::getState().
   * @return std::shared_ptr<Core::State> The current state
   */
  std::shared_ptr<Core::State> getState() const final;
  /**
   * @brief Implements Core::StateHandableObject::loadState().
   * @param state The new state.
   */
  void loadState(std::shared_ptr<Core::State> /*state*/) final;
  /**
   * @brief Whether the calculator has no underlying Python code and can therefore
   * release the global interpreter lock in Python bindings
   */
  bool allowsPythonGILRelease() const override {
    const auto underlyingCalculators = getUnderlyingCalculators();
    return std::all_of(underlyingCalculators.begin(), underlyingCalculators.end(),
                       [](const auto& c) { return c->allowsPythonGILRelease(); });
  };
  /**
   * @brief Getter for the MM calculator.
   * @return The MM calculator.
   */
  std::shared_ptr<MolecularMechanics::MolecularMechanicsCalculator> getMolecularMechanicsCalculator();
  /**
   * @brief Getter for the QM calculator.
   * @return The QM calculator.
   */
  std::shared_ptr<Core::Calculator> getQuantumMechanicsCalculator();

 private:
  void setStructureImpl(const Utils::AtomCollection& structure);
  void setUnderlyingCalculatorsImpl(std::vector<std::shared_ptr<Core::Calculator>> underlyingCalculators);
  void addUnderlyingSettingsImpl();
  void optimizeLinks();
  void handleElectrostaticEmbedding();
  void removeCalculatorSpecificSettings();
  /*
   * @brief Implementation of a calculation.
   */
  const Utils::Results& calculateImpl(std::string description);
  /*
   * @brief Implementation of a calculation where the QM calculation is not performed.
   */
  const Utils::Results& ignoreQmCalculateImpl(std::string description);
  /*
   * @brief Apply settings.
   */
  void applySettings();
  /*
   * @brief Apply settings to underlying calculators (QM and MM).
   */
  void applySettingsToUnderlyingCalculators();
  /*
   * @brief Eliminates those terms from the MM calculator which shall not be calculated at the MM level,
   *        but at the QM level instead.
   * @param reducedEnergyCalculation Whether the current calculation is part of a
   *                                 reduced QM/MM energy evaluation. Default: It is not.
   */
  void prepareTermsForMmCalculator(bool reducedEnergyCalculation = false);
  /*
   * @brief Sets the current log of this calculator to the underlying QM and MM calculators.
   */
  void setLogForUnderlyingCalculators();
  // The required properties.
  Utils::PropertyList requiredProperties_;
  // The settings.
  std::unique_ptr<Utils::Settings> settings_;
  // The results container.
  Utils::Results results_;
  // The structure of the full system
  Utils::AtomCollection structure_;
  // The structure of the QM region
  Utils::AtomCollection qmRegion_;
  // Vector containing the indices of the MM atoms close to the boundary (part of a cleft bond)
  std::vector<int> mmBoundaryAtoms_;
  // Vector containing the indices of the atoms in the QM region
  std::vector<int> listOfQmAtoms_;
  // QM calculator
  std::shared_ptr<Core::Calculator> qmCalculator_;
  // Flag for the case that the QM system spans the full system.
  bool mmAtomsLeft_ = true;
  /*
   * Whether an additional MM calculation is performed to evaluate the QM/MM energy
   * without any covalent and non-covalent contributions solely within the environment
   */
  bool calculateReducedQmMmEnergy_;
  // MM calculator
  std::shared_ptr<MolecularMechanics::MolecularMechanicsCalculator> mmCalculator_;
  // Whether electrostatic embedding shall be used (otherwise just mechanical embedding)
  bool electrostaticEmbedding_;
  // Whether to not perform the QM calculation.
  bool ignoreQm_;
  // File to which the QM region is written in XYZ format
  std::string qmRegionFile_;
  // The scheme for how to redistribute the charges close to the QM-MM boundary
  std::string chargeRedistributionScheme_;
  // Name of the point charges file
  static constexpr const char* pointChargesFilename_ = "environment_pointcharges.pc";
  /**
   * The following boolean keeps track of whether Turbomole is chosen as the QM calculator. If so, the formatting
   * of the point charges file requires a different file format.
   */
  bool turbomoleIsQmCalculator_ = false;
  /*
   * The following boolean keeps track of whether the SCF convergence setting for the QM calculator was already set,
   * because resetting it with every calculation (e.g., as part of an optimization) will result in a repeated
   * warning with some calculators (e.g., ORCA). The reason is that some calculators reset this setting internally
   * After calling setStructure(), this boolean will be reset to false again.
   *
   * Important note: This means that if one wants to reset the setting externally, setStructure() has to be recalled
   * so that the setting is passed on to the QM calculator.
   */
  bool scfConvCriterionIsSet_ = false;
  bool optimizeLinks_ = false;
};

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMMMCALCULATOR_H
