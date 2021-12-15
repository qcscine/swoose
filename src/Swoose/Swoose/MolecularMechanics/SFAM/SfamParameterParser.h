/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_SFAMPARAMETERPARSER_H
#define MOLECULARMECHANICS_SFAMPARAMETERPARSER_H

#include <istream>
#include <memory>
#include <regex>
#include <string>
#include <vector>

namespace Scine {
namespace MolecularMechanics {
class SfamParameters;
class AtomTypesHolder;

/**
 * @class SfamParameterParser SfamParameterParser.h
 * @brief This class reads in the SFAM parameters from the parameter file.
 *        The parsing relies on Qt and the conversion to doubles works also when numbers in the current locale are of
 *        the form "x,xxx".
 */
class SfamParameterParser {
 public:
  /**
   * @brief Constructor from filename of parameter file and atom types of the system.
   */
  SfamParameterParser(std::string filename, const AtomTypesHolder& atomTypes);

  /**
   * @brief Parse the parameters.
   */
  std::unique_ptr<SfamParameters> parseParameters();

 private:
  bool parse(SfamParameters& parameters);
  bool parseBonds(std::istream& in, SfamParameters& parameters);
  bool parseAngles(std::istream& in, SfamParameters& parameters);
  bool parseDihedrals(std::istream& in, SfamParameters& parameters);
  bool parseImproperDihedrals(std::istream& in, SfamParameters& parameters);
  bool parseCharges(std::istream& in, SfamParameters& parameters);
  bool parseNonCovalentParameters(std::istream& in, SfamParameters& parameters);
  bool parseC6Parameters(std::istream& in, SfamParameters& parameters);

  // Checks whether a regex iterator is valid; throws exception if not.
  void checkIter(const std::sregex_token_iterator& iter);

  std::string parameterFile_;
  int nAtoms_;
  const AtomTypesHolder& atomTypes_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_SFAMPARAMETERPARSER_H
