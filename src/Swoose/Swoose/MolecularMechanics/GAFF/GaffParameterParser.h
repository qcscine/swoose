/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_GAFFPARAMETERPARSER_H
#define MOLECULARMECHANICS_GAFFPARAMETERPARSER_H

#include <istream>
#include <memory>
#include <regex>
#include <string>
#include <unordered_map>
#include <vector>

namespace Scine {
namespace MolecularMechanics {
class GaffParameters;

/**
 * @class GaffParameterParser GaffParameterParser.h
 * @brief This class reads in the GAFF parameters from the parameter file.
 *        The parsing relies on Qt and the conversion to doubles works also when numbers in the current locale are of
 *        the form "x,xxx".
 */
class GaffParameterParser {
 public:
  /**
   * @brief Constructor from filename of parameter file.
   */
  explicit GaffParameterParser(std::string filename);

  /**
   * @brief Parse the parameters.
   */
  std::unique_ptr<GaffParameters> parseParameters();

 private:
  void parse(GaffParameters& parameters);
  void parseFirstLine(std::istream& in);
  void parseAtomTypesInfo(std::istream& in, GaffParameters& parameters);
  void parseHydrophilicAtomSymbols(std::istream& in, GaffParameters& parameters);
  void parseBonds(std::istream& in, GaffParameters& parameters);
  void parseAngles(std::istream& in, GaffParameters& parameters);
  void parseDihedrals(std::istream& in, GaffParameters& parameters);
  void parseImproperDihedrals(std::istream& in, GaffParameters& parameters);
  void parseHBond1012(std::istream& in, GaffParameters& parameters);
  void parseLennardJones(std::istream& in, GaffParameters& parameters);

  // Checks whether a regex iterator is valid; throws exception if not.
  void checkIter(const std::sregex_token_iterator& iter);

  std::string parameterFile_;
  std::unordered_map<std::string, double> atomicMasses_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_GAFFPARAMETERPARSER_H
