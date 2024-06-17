/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SfamParameterParser.h"
#include "../AtomTypesHolder.h"
#include "../Interactions/ElectrostaticEvaluator.h"
#include "../Interactions/RepulsionEvaluator.h"
#include "../Parameters/AngleParameters.h"
#include "../Parameters/BondParameters.h"
#include "SfamParameters.h"
#include <fstream>

namespace Scine {
namespace MolecularMechanics {

SfamParameterParser::SfamParameterParser(std::string filename, const AtomTypesHolder& atomTypes)
  : parameterFile_(std::move(filename)), nAtoms_(atomTypes.size()), atomTypes_(atomTypes) {
}

std::unique_ptr<SfamParameters> SfamParameterParser::parseParameters() {
  auto parameters = std::make_unique<SfamParameters>();

  bool successful = parse(*parameters);
  if (!successful)
    throw std::runtime_error("The parameter file is not valid!");

  return parameters;
}

bool SfamParameterParser::parse(SfamParameters& parameters) {
  std::ifstream file(parameterFile_);

  if (!file.is_open()) {
    throw std::runtime_error("The parameter file " + parameterFile_ + " cannot be opened.");
  }

  bool blockFound = parseBonds(file, parameters);
  if (!blockFound)
    return blockFound;

  blockFound = parseAngles(file, parameters);
  if (!blockFound)
    return blockFound;

  blockFound = parseDihedrals(file, parameters);
  if (!blockFound)
    return blockFound;

  blockFound = parseImproperDihedrals(file, parameters);
  if (!blockFound)
    return blockFound;

  blockFound = parseCharges(file, parameters);
  if (!blockFound)
    return blockFound;

  blockFound = parseNonCovalentParameters(file, parameters);
  if (!blockFound)
    return blockFound;

  blockFound = parseC6Parameters(file, parameters);
  return blockFound;
}

bool SfamParameterParser::parseBonds(std::istream& in, SfamParameters& parameters) {
  std::string keyword = "bond";
  std::string line;
  while (line.find(keyword) == std::string::npos) {
    if (!std::getline(in, line))
      return false;
  }
  if (!std::getline(in, line))
    return false;

  while (!line.empty() && line[0] != '!' && line[0] != '*') {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    if (*iter == "")
      iter++;

    checkIter(iter);
    std::string atomType1 = *iter++;
    checkIter(iter);
    std::string atomType2 = *iter++;

    checkIter(iter);
    double equilibriumDistance = std::stod(*iter++); // Angstrom
    checkIter(iter);
    double forceConstant = std::stod(*iter++); // kcal/mol/(A**2)

    parameters.addBond(BondType(atomType1, atomType2), BondParameters(forceConstant, equilibriumDistance));
    if (!std::getline(in, line))
      return false;
  }
  return true;
}

bool SfamParameterParser::parseAngles(std::istream& in, SfamParameters& parameters) {
  std::string keyword = "angle";
  std::string line;
  while (line.find(keyword) == std::string::npos) {
    if (!std::getline(in, line))
      return false;
  }
  if (!std::getline(in, line))
    return false;

  while (!line.empty() && line[0] != '!' && line[0] != '*') {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    if (*iter == "")
      iter++;

    checkIter(iter);
    std::string atomType1 = *iter++;
    checkIter(iter);
    std::string atomType2 = *iter++;
    checkIter(iter);
    std::string atomType3 = *iter++;

    checkIter(iter);
    double equilibriumAngle = std::stod(*iter++); // degrees
    checkIter(iter);
    double forceConstant = std::stod(*iter++); // kcal/mol / (rad^2)

    parameters.addAngle(AngleType(atomType1, atomType2, atomType3), AngleParameters(forceConstant, equilibriumAngle));
    if (!std::getline(in, line))
      return false;
  }
  return true;
}

bool SfamParameterParser::parseDihedrals(std::istream& in, SfamParameters& parameters) {
  std::string keyword = "dihedrals";
  std::string line;
  while (line.find(keyword) == std::string::npos) {
    if (!std::getline(in, line))
      return false;
  }
  if (!std::getline(in, line))
    return false;

  while (!line.empty() && line[0] != '!' && line[0] != '*') {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    if (*iter == "")
      iter++;

    checkIter(iter);
    std::string atomType1 = *iter++;
    checkIter(iter);
    std::string atomType2 = *iter++;
    checkIter(iter);
    std::string atomType3 = *iter++;
    checkIter(iter);
    std::string atomType4 = *iter++;

    checkIter(iter);
    double halfBarrier = std::stod(*iter++); // kcal/mol
    checkIter(iter);
    double phaseShift = std::stod(*iter++); // degrees
    checkIter(iter);
    int periodicity = std::stoi(*iter++);

    parameters.addDihedral(DihedralType(atomType1, atomType2, atomType3, atomType4),
                           DihedralParameters(halfBarrier, phaseShift, periodicity));
    if (!std::getline(in, line))
      return false;
  }
  return true;
}

bool SfamParameterParser::parseImproperDihedrals(std::istream& in, SfamParameters& parameters) {
  std::string keyword = "impropers";
  std::string line;
  while (line.find(keyword) == std::string::npos) {
    if (!std::getline(in, line))
      return false;
  }
  if (!std::getline(in, line))
    return false;

  while (!line.empty() && line[0] != '!' && line[0] != '*') {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    if (*iter == "")
      iter++;

    checkIter(iter);
    std::string atomType1 = *iter++; // central Atom
    checkIter(iter);
    std::string atomType2 = *iter++;
    checkIter(iter);
    std::string atomType3 = *iter++;
    checkIter(iter);
    std::string atomType4 = *iter++;

    checkIter(iter);
    double equilibriumAngle = std::stod(*iter++); // degrees
    checkIter(iter);
    double forceConstant = std::stod(*iter++); // kcal/mol

    parameters.addImproperDihedral(ImproperDihedralType(atomType1, atomType2, atomType3, atomType4),
                                   ImproperDihedralParameters(forceConstant, equilibriumAngle));
    if (!std::getline(in, line))
      return false;
  }
  return true;
}

bool SfamParameterParser::parseCharges(std::istream& in, SfamParameters& parameters) {
  std::string keyword = "charges";
  std::string line;
  while (line.find(keyword) == std::string::npos) {
    if (!std::getline(in, line))
      return false;
  }
  if (!std::getline(in, line))
    return false;

  while (!line.empty() && line[0] != '!' && line[0] != '*') {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    if (*iter == "")
      iter++;

    checkIter(iter);
    std::string atomTypeString = *iter++;
    checkIter(iter);
    double charge = std::stod(*iter++);

    parameters.addCharge(atomTypeString, charge);
    if (!std::getline(in, line))
      return false;
  }
  return true;
}

bool SfamParameterParser::parseNonCovalentParameters(std::istream& in, SfamParameters& parameters) {
  std::string keywordOne = "non-covalent";
  std::string keywordTwo = "vdw";
  std::string line;
  while (line.find(keywordOne) == std::string::npos && line.find(keywordTwo) == std::string::npos) {
    if (!std::getline(in, line))
      return false;
  }
  if (!std::getline(in, line))
    return false;

  std::vector<double> nonCovParams;
  while (!line.empty() && line[0] != '!' && line[0] != '*') {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);

    checkIter(iter);
    if (*iter == "")
      iter++;

    checkIter(iter);
    std::string parameterName = *iter++;
    checkIter(iter);
    double param = std::stod(*iter++);

    nonCovParams.push_back(param);
    if (!std::getline(in, line))
      return false;
  }

  // Add default values of the following two parameters if they were omitted in the parameter file
  if (nonCovParams.size() == 3)
    nonCovParams.push_back(Defaults::defaultBetaRepulsionParameter);
  if (nonCovParams.size() == 4)
    nonCovParams.push_back(Defaults::defaultAtomicChargesScalingFactor);

  parameters.setNonCovalentParameters(nonCovParams);
  return true;
}

bool SfamParameterParser::parseC6Parameters(std::istream& in, SfamParameters& parameters) {
  std::string keyword = "c6";
  std::string line;
  while (line.find(keyword) == std::string::npos) {
    if (!std::getline(in, line))
      return false;
  }
  if (!std::getline(in, line))
    return false;

  parameters.prepareC6Matrix(atomTypes_);
  int a = 0;
  while (!line.empty() && line[0] != '!' && line[0] != '*') {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
    std::sregex_token_iterator endOfLine;

    checkIter(iter);
    if (*iter == "")
      iter++;

    if (a >= int(parameters.getC6IndicesMap().size()))
      throw std::runtime_error("Error while parsing C6 coefficients from parameter file!");
    for (int b = 0; b <= a; ++b) {
      if (iter != endOfLine) {
        float c6 = std::stof(*iter++);
        parameters.setC6(a, b, c6);
      }
      else {
        throw std::runtime_error("Error while parsing C6 coefficients from parameter file!");
      }
    }
    if (iter != endOfLine) // no more values allowed after that in the current line
      throw std::runtime_error("Error while parsing C6 coefficients from parameter file!");
    a++;
    if (!std::getline(in, line))
      break; // at the very end of the file, no empty line or "*" symbol is needed!
  }
  // If there was not the correct number of rows parsed, reset the C6 matrix,
  // such that it can be detected from the outside (by the MM parameter sanity check)
  // that the parsing was not successful.
  if (a != int(parameters.getC6IndicesMap().size()))
    parameters.resetC6Matrix();
  return true;
}

void SfamParameterParser::checkIter(const std::sregex_token_iterator& iter) {
  std::sregex_token_iterator endOfLine;
  if (iter == endOfLine)
    throw std::runtime_error("The SFAM parameter file is not valid!");
}

} // namespace MolecularMechanics
} // namespace Scine
