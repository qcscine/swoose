/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "SwooseModule.h"
#include "MMParametrization/Parametrizer.h"
#include "MolecularMechanics/GAFF/GaffMolecularMechanicsCalculator.h"
#include "MolecularMechanics/SFAM/SfamMolecularMechanicsCalculator.h"
#include "QMMM/QmmmCalculator.h"
#include "StructurePreparation/StructureProcessor.h"
#include <Core/DerivedModule.h>
#include <Core/Exceptions.h>

namespace Scine {
namespace Swoose {

using InterfaceModelMap =
    boost::mpl::map<boost::mpl::pair<Core::Calculator, boost::mpl::vector<Qmmm::QmmmCalculator, MolecularMechanics::SfamMolecularMechanicsCalculator,
                                                                          MolecularMechanics::GaffMolecularMechanicsCalculator>>,
                    boost::mpl::pair<Core::MMParametrizer, boost::mpl::vector<MMParametrization::Parametrizer>>,
                    boost::mpl::pair<Core::EmbeddingCalculator, boost::mpl::vector<Qmmm::QmmmCalculator>>>;

std::string SwooseModule::name() const noexcept {
  return "Swoose";
}

boost::any SwooseModule::get(const std::string& interface, const std::string& model) const {
  boost::any resolved = Scine::Core::DerivedModule::resolve<InterfaceModelMap>(interface, model);

  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Scine::Core::ClassNotImplementedError();
  }

  return resolved;
}

bool SwooseModule::has(const std::string& interface, const std::string& model) const noexcept {
  return Scine::Core::DerivedModule::has<InterfaceModelMap>(interface, model);
}

std::vector<std::string> SwooseModule::announceInterfaces() const noexcept {
  return Scine::Core::DerivedModule::announceInterfaces<InterfaceModelMap>();
}

std::vector<std::string> SwooseModule::announceModels(const std::string& interface) const noexcept {
  return Scine::Core::DerivedModule::announceModels<InterfaceModelMap>(interface);
}

std::shared_ptr<Scine::Core::Module> SwooseModule::make() {
  return std::make_shared<SwooseModule>();
}

std::vector<std::shared_ptr<Scine::Core::Module>> moduleFactory() {
  return {SwooseModule::make()};
}

} // namespace Swoose
} // namespace Scine
