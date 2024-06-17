/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MolecularMechanicsCalculator.h"

namespace Scine {
namespace MolecularMechanics {

bool MolecularMechanicsCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == this->name();
}

std::unique_ptr<Utils::AtomCollection> MolecularMechanicsCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(structure_);
}

void MolecularMechanicsCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  structure_.setPositions(std::move(newPositions));
}

const Utils::PositionCollection& MolecularMechanicsCalculator::getPositions() const {
  return structure_.getPositions();
}

void MolecularMechanicsCalculator::setRequiredProperties(const Utils::PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
}

Utils::PropertyList MolecularMechanicsCalculator::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian | Utils::Property::AtomicCharges |
         Utils::Property::SuccessfulCalculation | Utils::Property::BondOrderMatrix | Utils::Property::PartialEnergies;
}

Utils::PropertyList MolecularMechanicsCalculator::getRequiredProperties() const {
  return requiredProperties_;
}

const Utils::Settings& MolecularMechanicsCalculator::settings() const {
  return *settings_;
}

Utils::Settings& MolecularMechanicsCalculator::settings() {
  return *settings_;
}

const Utils::Results& MolecularMechanicsCalculator::results() const {
  return results_;
}

Utils::Results& MolecularMechanicsCalculator::results() {
  return results_;
}

const std::vector<double>& MolecularMechanicsCalculator::atomicCharges() const {
  return atomicCharges_;
}

const std::vector<std::list<int>>& MolecularMechanicsCalculator::listsOfNeighbors() const {
  return listsOfNeighbors_;
}

void MolecularMechanicsCalculator::loadState(std::shared_ptr<Core::State> /*state*/) {
  // empty implementation
}

std::shared_ptr<Core::State> MolecularMechanicsCalculator::getState() const {
  return nullptr;
}

} // namespace MolecularMechanics
} // namespace Scine
