/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CommandLineOptions.h"
#include <Core/Exceptions.h>
#include <Core/Log.h>
#include <boost/program_options.hpp>

namespace Scine {
namespace Swoose {

namespace {
constexpr const char* helpKey = "help";
constexpr const char* settingsFileKey = "settings";
constexpr const char* structureKey = "structure";
constexpr const char* modeKey = "mode";
constexpr const char* calculateHessianKey = "hessian";
constexpr const char* logDebugOptionKey = "verbose";
constexpr const char* quantumKey = "quantum";
} // namespace

using namespace boost::program_options;

struct CommandLineOptions::Impl {
  boost::program_options::options_description desc_;
  boost::program_options::variables_map vm_;
  std::string callStatement;
};

CommandLineOptions::CommandLineOptions(int argc, char* argv[]) : pImpl_(std::make_unique<Impl>()) {
  pImpl_->callStatement = generateCallStatement(argc, argv);

  auto helpOption = combineNamesForOptions(helpKey, ",h");
  auto modeOption = combineNamesForOptions(modeKey, ",m");
  auto settingsFileOption = combineNamesForOptions(settingsFileKey, ",y");
  auto structureOption = combineNamesForOptions(structureKey, ",s");
  auto quantumOption = combineNamesForOptions(quantumKey, ",q");
  auto hessianOption = combineNamesForOptions(calculateHessianKey, ",H");
  auto logDebugOption = combineNamesForOptions(logDebugOptionKey, ",v");

  // clang-format off
  pImpl_->desc_.add_options()
    (helpOption.c_str(), "Prints this help message.")
    (modeOption.c_str(), value<std::string>(), "Sets the mode of the Swoose app. Options: parametrize, calculate, md, optimize, select_qm.")
    (settingsFileOption.c_str(), value<std::string>(), "Sets the path to the settings YAML file.")
    (structureOption.c_str(), value<std::string>(), "Sets the path to the molecular structure's XYZ file.")
    (quantumOption.c_str(), "Activates the QM/MM hybrid model.")
    (hessianOption.c_str(), "Activates the calculation of the Hessian.")
    (logDebugOption.c_str(), "Activates a more verbose output.");

  // clang-format on

  positional_options_description p;
  store(command_line_parser(argc, argv).options(pImpl_->desc_).positional(p).run(), pImpl_->vm_);
  notify(pImpl_->vm_);
}

CommandLineOptions::~CommandLineOptions() = default;

std::string CommandLineOptions::generateCallStatement(int argc, char* argv[]) const {
  std::string callStatement = "";
  for (int i = 0; i < argc; ++i) {
    callStatement += " " + std::string(argv[i]);
  }
  return callStatement;
}

template<class CharPtrType, class StringType>
std::string CommandLineOptions::combineNamesForOptions(CharPtrType nameOfOption, StringType abbreviatedOption) const {
  auto combinedString = static_cast<std::string>(nameOfOption) + static_cast<std::string>(abbreviatedOption);
  return combinedString;
}

bool CommandLineOptions::helpRequired() const {
  return pImpl_->vm_.count(helpKey) > 0;
}

void CommandLineOptions::printHelp(Core::Log& log) const {
  log.output << "Usage: <executable> [options]\n";
  log.output << pImpl_->desc_;
}

std::string CommandLineOptions::getStructureFile() const {
  if (pImpl_->vm_.count(structureKey) > 0) {
    return pImpl_->vm_[structureKey].as<std::string>();
  }
  throw Core::InitializationException("No molecular structure file given.");
}

std::string CommandLineOptions::getSettingsFile() const {
  if (pImpl_->vm_.count(settingsFileKey) > 0) {
    return pImpl_->vm_[settingsFileKey].as<std::string>();
  }
  return "";
}

std::string CommandLineOptions::getMode() const {
  if (pImpl_->vm_.count(modeKey) > 0) {
    return pImpl_->vm_[modeKey].as<std::string>();
  }
  throw Core::InitializationException("No mode was set.");
}

bool CommandLineOptions::hessianRequired() const {
  return pImpl_->vm_.count(calculateHessianKey) > 0;
}

bool CommandLineOptions::quantumCalculationRequired() const {
  return pImpl_->vm_.count(quantumKey) > 0;
}

bool CommandLineOptions::debugLoggingRequired() const {
  return pImpl_->vm_.count(logDebugOptionKey) > 0;
}

} // namespace Swoose
} // namespace Scine
