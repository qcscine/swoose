/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaffParameterParser.h"
#include "../AtomTypesHolder.h"
#include "../Parameters/AngleParameters.h"
#include "../Parameters/BondParameters.h"
#include "../Parameters/LennardJonesParameters.h"
#include "GaffParameters.h"
#include <fstream>

namespace Scine {
namespace MolecularMechanics {

GaffParameterParser::GaffParameterParser(std::string filename) : parameterFile_(std::move(filename)) {
}

std::unique_ptr<GaffParameters> GaffParameterParser::parseParameters() {
  auto parameters = std::make_unique<GaffParameters>();
  parse(*parameters);
  return parameters;
}

void GaffParameterParser::parse(GaffParameters& parameters) {
  std::ifstream file(parameterFile_);

  if (!file.is_open()) {
    throw std::runtime_error("The parameter file " + parameterFile_ + " cannot be opened.");
  }

  parseFirstLine(file);
  parseAtomTypesInfo(file, parameters);
  parseHydrophilicAtomSymbols(file, parameters);
  parseBonds(file, parameters);
  parseAngles(file, parameters);
  parseDihedrals(file, parameters);
  parseImproperDihedrals(file, parameters);
  parseHBond1012(file, parameters);
  parseLennardJones(file, parameters);
}

void GaffParameterParser::parseFirstLine(std::istream& in) {
  // NB: first line contains general info about parameters, it is not needed
  std::string line;
  std::getline(in, line);
}

void GaffParameterParser::parseAtomTypesInfo(std::istream& in, GaffParameters& /*parameters*/) {
  std::string line;
  std::getline(in, line);
  while (!line.empty()) {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    std::string atomType = *iter++;
    checkIter(iter);
    std::string massString = *iter++;
    checkIter(iter);
    std::string atomicPolarizability = *iter++; // NB: not needed

    // atomicMasses_[atomType] = std::stod(massString); // NB: also not needed
    std::getline(in, line);
  }
}

void GaffParameterParser::parseHydrophilicAtomSymbols(std::istream& in, GaffParameters& /*parameters*/) {
  // NB: line containing list of hydrophylic atoms. Not needed
  std::string line;
  std::getline(in, line);
}

void GaffParameterParser::parseBonds(std::istream& in, GaffParameters& parameters) {
  std::string line;
  std::getline(in, line);
  while (!line.empty()) {
    std::regex rgx("(\\s|-)+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    std::string atomType1 = *iter++;
    checkIter(iter);
    std::string atomType2 = *iter++;
    checkIter(iter);
    std::string forceConstantString = *iter++;
    checkIter(iter);
    std::string equilibriumDistanceString = *iter++;
    double forceConstant = std::stod(forceConstantString);             // kcal/mol/(A**2)
    double equilibriumDistance = std::stod(equilibriumDistanceString); // Angstrom

    forceConstant *= 2.0; // In GAFF, the harmonic potential is k*x^2 instead of 0.5*k*x^2.

    parameters.addBond(BondType(atomType1, atomType2), BondParameters(forceConstant, equilibriumDistance));
    std::getline(in, line);
  }
}

void GaffParameterParser::parseAngles(std::istream& in, GaffParameters& parameters) {
  std::string line;
  std::getline(in, line);
  while (!line.empty()) {
    std::regex rgx("(\\s|-)+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    std::string atomType1 = *iter++;
    checkIter(iter);
    std::string atomType2 = *iter++;
    checkIter(iter);
    std::string atomType3 = *iter++;
    checkIter(iter);
    double forceConstant = std::stod(*iter++); // kcal/mol / (rad^2)
    checkIter(iter);
    double equilibriumAngle = std::stod(*iter++); // degrees

    forceConstant *= 2.0; // In GAFF, the harmonic potential is k*x^2 instead of 0.5*k*x^2.

    parameters.addAngle(AngleType(atomType1, atomType2, atomType3), AngleParameters(forceConstant, equilibriumAngle));
    std::getline(in, line);
  }
}

void GaffParameterParser::parseDihedrals(std::istream& in, GaffParameters& parameters) {
  std::string line;
  std::getline(in, line);
  while (!line.empty()) {
    std::regex rgx("(\\s|-)+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    std::string atomType1 = *iter++;
    checkIter(iter);
    std::string atomType2 = *iter++;
    checkIter(iter);
    std::string atomType3 = *iter++;
    checkIter(iter);
    std::string atomType4 = *iter++;

    checkIter(iter);
    unsigned divisionFactor = static_cast<unsigned>(std::stoi(*iter++));
    checkIter(iter);
    double halfBarrier = std::stod(*iter++); // kcal/mol
    checkIter(iter);
    double phaseShift = std::stod(*iter++); // degrees
    checkIter(iter);
    // NB: minus sign in case of multiple dihedrals are removed by regex (see above)
    unsigned periodicity = static_cast<unsigned>(std::lround(std::stod(*iter++)));

    halfBarrier /= divisionFactor;
    parameters.addDihedral(DihedralType(atomType1, atomType2, atomType3, atomType4),
                           DihedralParameters(halfBarrier, phaseShift, periodicity));
    std::getline(in, line);
  }
}

void GaffParameterParser::parseImproperDihedrals(std::istream& in, GaffParameters& parameters) {
  std::string line;
  std::getline(in, line);
  while (!line.empty()) {
    std::regex rgx("(\\s|-)+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    std::string atomType1 = *iter++;
    checkIter(iter);
    std::string atomType2 = *iter++;
    checkIter(iter);
    std::string atomType3 = *iter++;
    checkIter(iter);
    std::string atomType4 = *iter++;

    // unsigned divisionFactor = static_cast<unsigned>(std::stoi(*iter++)); // not specified for improper dihedrals
    checkIter(iter);
    double halfBarrier = std::stod(*iter++); // kcal/mol
    checkIter(iter);
    double phaseShift = std::stod(*iter++); // degrees
    checkIter(iter);
    // NB: minus sign in case of multiple dihedrals are removed by regex (see above)
    unsigned periodicity = static_cast<unsigned>(std::lround(std::stod(*iter++)));

    // In GAFF, the improper dihedrals are essentially treated as normal dihedrals
    parameters.addImproperDihedral(ImproperDihedralType(atomType3, atomType1, atomType2, atomType4),
                                   DihedralParameters(halfBarrier, phaseShift, periodicity));
    std::getline(in, line);
  }
}

void GaffParameterParser::parseHBond1012(std::istream& in, GaffParameters& /*parameters*/) {
  // NB: not handled yet (is zero in gaff)
  std::string line;
  std::getline(in, line);
  while (!line.empty()) {
    std::getline(in, line);
  }
}

void GaffParameterParser::parseLennardJones(std::istream& in, GaffParameters& parameters) {
  std::string line;
  std::getline(in, line);

  // First look at which atoms are to be equivalenced
  std::vector<std::pair<std::string, std::string>> equivalentLennardJonesAtoms;
  while (!line.empty()) {
    std::regex rgx("(\\s|-)+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    std::string firstAtom = *iter++;
    std::sregex_token_iterator end;
    for (; iter != end; ++iter) {
      checkIter(iter);
      equivalentLennardJonesAtoms.emplace_back(*iter, firstAtom);
    }

    std::getline(in, line);
  }

  // Reads label and kind
  // Label is irrelevant for us, as I understand the main program loads the parameters by looking for labels in this
  // file. Kind: so far only RE is supported
  std::getline(in, line);
  std::regex rgx("\\s+");
  std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
  checkIter(iter);
  std::string nonBondedIgnoredString = *iter++;
  checkIter(iter);
  std::string re = *iter++;

  if (re != "RE")
    throw std::runtime_error("GAFF non-bonded format not supported.");

  // Read non-bonded parameters ("RE")
  std::getline(in, line);
  while (!line.empty()) {
    std::regex rgx2("\\s+");
    std::sregex_token_iterator iter2(line.begin(), line.end(), rgx2, -1);

    checkIter(iter2);
    *iter2++; // ignore preceding spaces
    checkIter(iter2);
    std::string atom = *iter2++;
    checkIter(iter2);
    double vdwRadius = std::stod(*iter2++); // Angstrom
    checkIter(iter2);
    double wellDepth = std::stod(*iter2++); // kcal/mol

    parameters.addLennardJones(atom, LennardJonesParameters(vdwRadius, wellDepth));

    // Look if other atom types were declared to have the same parameters:
    for (const auto& p : equivalentLennardJonesAtoms) {
      if (p.second == atom)
        parameters.addLennardJones(p.first, LennardJonesParameters(vdwRadius, wellDepth));
    }
    std::getline(in, line);
  }
}

void GaffParameterParser::checkIter(const std::sregex_token_iterator& iter) {
  std::sregex_token_iterator endOfLine;
  if (iter == endOfLine)
    throw std::runtime_error("The GAFF parameter file is not valid!");
}

} // namespace MolecularMechanics
} // namespace Scine
