/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_MOCK_QMCALCULATOR_H
#define SWOOSE_MOCK_QMCALCULATOR_H

#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {
namespace Swoose {

namespace SettingsNames {
static constexpr const char* mockIntSettingKey = "mock_int_setting";
static constexpr const char* mockStringSettingKey = "mock_string_setting";
} // namespace SettingsNames

/// @brief Mocks the settings for the mock QM calculator.
class MockQmCalculatorSettings : public Scine::Utils::Settings {
 public:
  MockQmCalculatorSettings() : Utils::Settings("MockQmCalculatorSettings") {
    Utils::UniversalSettings::IntDescriptor mockIntSetting("Mock integer setting.");
    mockIntSetting.setDefaultValue(42);
    _fields.push_back(SettingsNames::mockIntSettingKey, std::move(mockIntSetting));

    Utils::UniversalSettings::StringDescriptor mockStringSetting("Mock string setting.");
    mockStringSetting.setDefaultValue("test");
    _fields.push_back(SettingsNames::mockStringSettingKey, std::move(mockStringSetting));

    Utils::UniversalSettings::IntDescriptor spinMultiplicity("Spin multiplicity.");
    spinMultiplicity.setDefaultValue(1);
    _fields.push_back(Utils::SettingsNames::spinMultiplicity, std::move(spinMultiplicity));

    Utils::UniversalSettings::IntDescriptor molecularCharge("Molecular charge.");
    molecularCharge.setDefaultValue(0);
    _fields.push_back(Utils::SettingsNames::molecularCharge, std::move(molecularCharge));

    resetToDefaults();
  }
};

/// @brief Mock class for a QM calculator to use in QM/MM tests.
class MockQmCalculator : public Utils::CloneInterface<MockQmCalculator, Core::Calculator> {
 public:
  static constexpr const char* model = "MOCK-QM";
  /// @brief Default Constructor.
  MockQmCalculator();
  /// @brief Default Destructor.
  ~MockQmCalculator() final = default;
  /// @brief Copy Constructor.
  MockQmCalculator(const MockQmCalculator& rhs);
  /**
   * @brief Changes the molecular structure to calculate.
   * @param structure A new AtomCollection to save.
   */
  void setStructure(const Utils::AtomCollection& structure) override;
  /**
   * @brief Gets the molecular structure as a std::unique_ptr<AtomCollection>.
   * @return std::unique_ptr<AtomCollection>
   */
  std::unique_ptr<Utils::AtomCollection> getStructure() const override;
  /**
   * @brief Allows to modify the positions of the underlying AtomCollection
   * @param newPositions the new positions to be assigned to the underlying AtomCollection
   */
  void modifyPositions(Utils::PositionCollection newPositions) override;
  /**
   * @brief Getter for the coordinates of the underlying AtomCollection
   */
  const Utils::PositionCollection& getPositions() const override;
  /**
   * @brief Sets the properties to calculate.
   * @param requiredProperties A PropertyList, a sequence of bits that represent the
   *        properties that must be calculated.
   */
  void setRequiredProperties(const Utils::PropertyList& requiredProperties) override;
  /**
   * @brief Getter for the properties to calculate.
   */
  Utils::PropertyList getRequiredProperties() const override;
  /**
   * @brief Returns the list of the possible properties to calculate analytically.
   * By some method analytical hessian calculation is not possible. In this case the
   * hessian calculation is done seminumerically.
   */
  Utils::PropertyList possibleProperties() const override;
  /**
   * @brief The main function running calculations (dummy).
   *
   * @param description   The calculation description.
   * @return Result Return the result of the calculation. The object contains the
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
   * @return Settings& The settings.
   */
  Utils::Settings& settings() override;
  /**
   * @brief Constant accessor for the settings.
   * @return const Settings& The settings.
   */
  const Utils::Settings& settings() const override;
  /**
   * @brief Implements Core::StateHandableObject::getState().
   * @return std::shared_ptr<Core::State> The current state
   */
  std::shared_ptr<Core::State> getState() const final;
  /**
   * @brief Implements Core::StateHandableObject::loadState().
   * @param state The new state.
   */
  void loadState(std::shared_ptr<Core::State> state) final;
  /**
   * @brief Accessor for the saved instance of Results.
   * @return Results& The results of the previous calculation.
   */
  Utils::Results& results() override;
  /**
   * @brief Constant accessor for the Results.
   * @return const Results& The results of the previous calculation.
   */
  const Utils::Results& results() const override;
  /**
   * @brief Getter for the file name base string.
   */
  std::string getFileNameBase() const;
  /**
   * @brief Getter for the calculation directory.
   */
  std::string getCalculationDirectory() const;
  /**
   * @brief Whether the calculator supports a method family
   * @param methodFamily identifier for the method family
   * @return whether the calculator supports a method family
   */
  bool supportsMethodFamily(const std::string& methodFamily) const override;

 private:
  /*
   * @brief Apply settings.
   */
  void applySettings();
  // The settings.
  std::unique_ptr<Utils::Settings> settings_;
  // The results.
  Utils::Results results_;
  // The required properties.
  Utils::PropertyList requiredProperties_;
  // The molecular structure.
  Utils::AtomCollection structure_;
  // An integer setting.
  int intSetting_;
  // A string setting.
  std::string stringSetting_;
};

} // namespace Swoose
} // namespace Scine

#endif // SWOOSE_MOCK_QMCALCULATOR_H
