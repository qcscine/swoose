/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_COMMANDLINEOPTIONS_H
#define SWOOSE_COMMANDLINEOPTIONS_H

#include <memory>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace Swoose {

/**
 *
 * @class CommandLineOptions CommandLineOptions.h
 * @brief Class to parse the command line options for non-default options and pass them to a Utils::Settings class.
 *
 *        This class uses the pImpl idiom to hide the boost::program_options dependency.
 */
class CommandLineOptions {
 public:
  /**
   * @brief Class constructor, parses the command line arguments and maps them to the according setting.
   * @param argv the vector with the argument strings.
   * @param argc the number of arguments.
   */
  CommandLineOptions(int argc, char* argv[]);
  ~CommandLineOptions();

  /** @brief Returns whether the help flag option has been set. */
  bool helpRequired() const;
  /** @brief Prints the help message. */
  void printHelp(Core::Log& log) const;

  /** @brief Returns the path to the molecular structure's XYZ file. */
  std::string getStructureFile() const;

  /** @brief Returns the path to the settings YAML file. */
  std::string getSettingsFile() const;

  /** @brief Returns the mode of the Swoose app. */
  std::string getMode() const;

  /** @brief Returns whether a QM/MM hybrid model is turned on. */
  bool quantumCalculationRequired() const;

  /** @brief Returns whether a Hessian is required from an MM calculation */
  bool hessianRequired() const;

  /** @brief Returns whether debug information should be logged. */
  bool debugLoggingRequired() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> pImpl_;

  /// Parses the command line to generate a call statement.
  std::string generateCallStatement(int argc, char* argv[]) const;

  /// Combines two character types to form a single const char*.
  template<class CharPtrType, class StringType>
  std::string combineNamesForOptions(CharPtrType nameOfOption, StringType abbreviatedOption) const;
};

} // namespace Swoose
} // namespace Scine

#endif // SWOOSE_COMMANDLINEOPTIONS_H
