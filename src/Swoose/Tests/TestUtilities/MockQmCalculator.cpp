/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MockQmCalculator.h"

namespace Scine {
namespace Swoose {

std::string MockQmCalculator::name() const {
  return "MOCK-QM";
}

bool MockQmCalculator::supportsMethodFamily(const std::string& /* methodFamily */) const {
  return true;
}

MockQmCalculator::MockQmCalculator() {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<MockQmCalculatorSettings>();
  applySettings();
}

MockQmCalculator::MockQmCalculator(const MockQmCalculator& rhs) {
  this->requiredProperties_ = rhs.requiredProperties_;
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  applySettings();
  this->results() = rhs.results();
  this->setStructure(rhs.structure_);
}

void MockQmCalculator::applySettings() {
  settings_->normalizeStringCases(); // convert all option names to lower case letters
  if (!settings_->valid()) {
    settings_->throwIncorrectSettings();
  }
  intSetting_ = settings_->getInt(SettingsNames::mockIntSettingKey);
  stringSetting_ = settings_->getString(SettingsNames::mockStringSettingKey);
}

void MockQmCalculator::setStructure(const Utils::AtomCollection& structure) {
  applySettings();
  structure_ = structure;
}

std::unique_ptr<Utils::AtomCollection> MockQmCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(structure_);
}

void MockQmCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  structure_.setPositions(std::move(newPositions));
}

const Utils::PositionCollection& MockQmCalculator::getPositions() const {
  return structure_.getPositions();
}

void MockQmCalculator::setRequiredProperties(const Utils::PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
}

Utils::PropertyList MockQmCalculator::getRequiredProperties() const {
  return requiredProperties_;
}

Utils::PropertyList MockQmCalculator::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients;
}

const Utils::Results& MockQmCalculator::calculate(std::string description) {
  applySettings();

  double energy = 200.0;
  Utils::GradientCollection gradients(structure_.size(), 3);
  for (int i = 0; i < structure_.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      gradients(i, j) = 100.0;
    }
  }

  results_.set<Utils::Property::Description>(description);
  if (requiredProperties_.containsSubSet(Utils::Property::Energy))
    results_.set<Utils::Property::Energy>(energy);
  if (requiredProperties_.containsSubSet(Utils::Property::Gradients))
    results_.set<Utils::Property::Gradients>(gradients);

  results_.set<Utils::Property::SuccessfulCalculation>(true);
  results_.set<Utils::Property::ProgramName>("mock");
  return results_;
}

const Utils::Settings& MockQmCalculator::settings() const {
  return *settings_;
}

Utils::Settings& MockQmCalculator::settings() {
  return *settings_;
}

std::shared_ptr<Core::State> MockQmCalculator::getState() const {
  return nullptr;
}

void MockQmCalculator::loadState(std::shared_ptr<Core::State> /* state */) {
  // Empty implementation
}

const Utils::Results& MockQmCalculator::results() const {
  return results_;
}

Utils::Results& MockQmCalculator::results() {
  return results_;
}

} // namespace Swoose
} // namespace Scine