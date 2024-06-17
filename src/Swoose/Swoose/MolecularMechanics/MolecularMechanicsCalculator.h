/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_MOLECULARMECHANICSCALCULATOR_H
#define MOLECULARMECHANICS_MOLECULARMECHANICSCALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <list>

namespace Scine {
namespace MolecularMechanics {

/**
 * @class MolecularMechanicsCalculator MolecularMechanicsCalculator.h
 * @brief Base class for the MM methods (currently: SFAM and GAFF).
 *
 */
class MolecularMechanicsCalculator
  : public Utils::CloneInterface<Utils::Abstract<MolecularMechanicsCalculator>, Core::Calculator> {
 public:
  /// @brief Default Constructor.
  MolecularMechanicsCalculator() = default;
  /// @brief Default Destructor.
  virtual ~MolecularMechanicsCalculator() override = default;
  MolecularMechanicsCalculator(const MolecularMechanicsCalculator& rhs) : CloneInterface(rhs){};
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
   * @brief Const accessor for the atomic charges.
   */
  const std::vector<double>& atomicCharges() const;
  /**
   * @brief Const accessor for the lists of neighbors (the connectivity of the molecular system).
   */
  const std::vector<std::list<int>>& listsOfNeighbors() const;
  /**
   * @brief Whether the calculator has no underlying Python code and can therefore
   * release the global interpreter lock in Python bindings
   */
  bool allowsPythonGILRelease() const override {
    return true;
  };

 protected:
  std::unique_ptr<Utils::Settings> settings_;
  Utils::AtomCollection structure_;
  Utils::PropertyList requiredProperties_;
  Utils::Results results_;
  std::vector<std::list<int>> listsOfNeighbors_;
  std::vector<double> atomicCharges_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_MOLECULARMECHANICSCALCULATOR_H
