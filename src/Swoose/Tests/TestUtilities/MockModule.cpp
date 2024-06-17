/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MockModule.h"
#include "MockQmCalculator.h"
#include <Core/DerivedModule.h>
#include <Swoose/MolecularMechanics/SFAM/SfamMolecularMechanicsCalculator.h>

namespace Scine {
namespace Swoose {

/*
 * Note: The SFAM calculator is added to this module such that a MM calculator
 *       is available in the test through the module system.
 */
using InterfaceModelMap =
    boost::mpl::map<boost::mpl::pair<Scine::Core::Calculator, boost::mpl::vector<MockQmCalculator, MolecularMechanics::SfamMolecularMechanicsCalculator>>>;

std::string MockModule::name() const noexcept {
  return "MockModule";
}

boost::any MockModule::get(const std::string& interface, const std::string& model) const {
  boost::any resolved = Scine::Core::DerivedModule::resolve<InterfaceModelMap>(interface, model);

  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }

  return resolved;
}

bool MockModule::has(const std::string& interface, const std::string& model) const noexcept {
  return Scine::Core::DerivedModule::has<InterfaceModelMap>(interface, model);
}

std::vector<std::string> MockModule::announceInterfaces() const noexcept {
  return Scine::Core::DerivedModule::announceInterfaces<InterfaceModelMap>();
}

std::vector<std::string> MockModule::announceModels(const std::string& interface) const noexcept {
  return Scine::Core::DerivedModule::announceModels<InterfaceModelMap>(interface);
}

std::shared_ptr<Scine::Core::Module> MockModule::make() {
  return std::make_shared<MockModule>();
}

std::vector<std::shared_ptr<Scine::Core::Module>> moduleFactory() {
  return {MockModule::make()};
}

} // namespace Swoose
} // namespace Scine
