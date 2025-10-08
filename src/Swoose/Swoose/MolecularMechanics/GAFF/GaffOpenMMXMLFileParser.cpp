/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Swoose/MolecularMechanics/GAFF/GaffOpenMMXMLFileParser.h"      // Class header file.
#include "Swoose/MolecularMechanics/GAFF/GaffParameters.h"               // Final result.
#include "Swoose/MolecularMechanics/Parameters/AngleParameters.h"        // Add angle parameters.
#include "Swoose/MolecularMechanics/Parameters/BondParameters.h"         // Add bond parameters.
#include "Swoose/MolecularMechanics/Parameters/LennardJonesParameters.h" // Add LJ parameters.
#include <Utils/Constants.h>                                             // Conversion between units.

namespace Scine {
namespace MolecularMechanics {

constexpr double angstrom_per_nanometer = 10.0;
constexpr double kcalPerMol_per_kJPerMol = Utils::Constants::kCalPerMol_per_hartree * Utils::Constants::hartree_per_kJPerMol;
constexpr double kcalPerMolPerAngstrom2_per_kJPerMolPerNanometer2 =
    kcalPerMol_per_kJPerMol / angstrom_per_nanometer / angstrom_per_nanometer;

GaffParameters GaffOpenMMXMLFileParser::parseParameters(const std::vector<std::string>& xmlFilePaths, Core::Log& log) {
  GaffParameters params;
  for (const auto& path : xmlFilePaths) {
    params += parseParameters(path, log);
    std::cout << "----------" << std::endl;
  }
  return params;
}

GaffParameters GaffOpenMMXMLFileParser::parseParameters(const std::string& xmlFilePath, Core::Log& log) {
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(xmlFilePath.c_str());
  GaffParameters parameters;
  if (result) {
    const auto type2class = parseType2AtomClassMap(doc);
    parameters += parseHarmonicBondForce(doc, log);
    parameters += parseHarmonicAngleForce(doc, log);
    parameters += parsePeriodicTorsionForce(doc, log);
    parameters += parseImproperPeriodicTorsionForce(doc, log);
    parameters += parseNonbondedForce(doc, log);
    parameters.setType2AtomClass(type2class);
  }
  else {
    throw std::runtime_error("Failed to load XML file:" + xmlFilePath +
                             ". PugiXML returned the following message: " + result.description());
  }
  if (parameters.empty()) {
    throw std::runtime_error("No parameters could be read from the XML file. The file may be empty or use an"
                             "unexpected notation.");
  }
  return parameters;
}

GaffParameters GaffOpenMMXMLFileParser::parseHarmonicBondForce(const pugi::xml_document& xmlDocument, Core::Log& log) {
  GaffParameters parameters;
  pugi::xml_node forces = xmlDocument.child("ForceField").child("HarmonicBondForce");

  for (pugi::xml_node force = forces.child("Bond"); force; force = force.next_sibling("Bond")) {
    std::string type1;
    std::string type2;
    std::string class1;
    std::string class2;
    double forceConstant = NAN;
    double length = NAN;
    try {
      type1 = force.attribute("type1").as_string();
      type2 = force.attribute("type2").as_string();
      class1 = force.attribute("class1").as_string();
      class2 = force.attribute("class2").as_string();
      if (force.attribute("k")) {
        forceConstant = force.attribute("k").as_double() * kcalPerMolPerAngstrom2_per_kJPerMolPerNanometer2; // kcal/mol/angstrom^2
      }
      if (force.attribute("length")) {
        length = force.attribute("length").as_double() * angstrom_per_nanometer; // angstrom
      }
    }
    catch (...) {
      throw std::runtime_error("Failed to read harmonic bond force definition in XML file.");
    }
    if (std::isnan(forceConstant) || std::isnan(length)) {
      throw std::runtime_error(std::string("Harmonic bond force information in XML file incomplete. The following was "
                                           "read from the file:\n") +
                               " type1=" + type1 + " , type2=" + type2 + " class1=" + class1 + " , class2=" + class2 +
                               " | length=" + std::to_string(length) + " | k=" + std::to_string(forceConstant));
    }
    if (!class1.empty() && !class2.empty()) {
      parameters.addBond(BondType(class1, class2), BondParameters(forceConstant, length));
    }
    else {
      if (type1.empty() || type2.empty()) {
        throw std::runtime_error("Harmonic bond information: Missing atom types!");
      }
      parameters.addBond(BondType(type1, type2), BondParameters(forceConstant, length));
    }
  }
  if (parameters.empty()) {
    log.warning << "Warning: "
                << "No HarmonicBondForces were read from the XML file!" << Core::Log::endl;
  }
  return parameters;
}

GaffParameters GaffOpenMMXMLFileParser::parseHarmonicAngleForce(const pugi::xml_document& xmlDocument, Core::Log& log) {
  GaffParameters parameters;
  pugi::xml_node forces = xmlDocument.child("ForceField").child("HarmonicAngleForce");
  for (pugi::xml_node force = forces.child("Angle"); force; force = force.next_sibling("Angle")) {
    std::string type1, class1;
    std::string type2, class2;
    std::string type3, class3;
    double forceConstant = NAN;
    double angle = NAN;
    try {
      type1 = force.attribute("type1").as_string();
      type2 = force.attribute("type2").as_string();
      type3 = force.attribute("type3").as_string();
      class1 = force.attribute("class1").as_string();
      class2 = force.attribute("class2").as_string();
      class3 = force.attribute("class3").as_string();
      if (force.attribute("k")) {
        forceConstant = force.attribute("k").as_double() * kcalPerMol_per_kJPerMol; // kcal/mol/rad^2
      }
      if (force.attribute("angle")) {
        angle = force.attribute("angle").as_double() * Utils::Constants::degree_per_rad; // degree
      }
    }
    catch (...) {
      throw std::runtime_error("Failed to read harmonic angle force definition in XML file.");
    }
    if (std::isnan(forceConstant) || std::isnan(angle)) {
      throw std::runtime_error(std::string("Harmonic angle force information in XML file incomplete. The following was "
                                           "read from the file:\n") +
                               " type1=" + type1 + " , type2=" + type2 + " , type3=" + type3 + " class1=" + class1 +
                               " , class2=" + class2 + " , class3=" + class3 + " | angle=" + std::to_string(angle) +
                               " | k=" + std::to_string(forceConstant));
    }
    if (!class1.empty() && !class2.empty() && !class3.empty()) {
      parameters.addAngle(AngleType(class1, class2, class3), AngleParameters(forceConstant, angle));
    }
    else {
      if (type1.empty() || type2.empty() || type3.empty()) {
        throw std::runtime_error("Harmonic angle information: Missing atom types!");
      }
      parameters.addAngle(AngleType(type1, type2, type3), AngleParameters(forceConstant, angle));
    }
  }
  if (parameters.empty()) {
    log.warning << "Warning: "
                << "No HarmonicAngleForces were read from the XML file!" << Core::Log::endl;
  }
  return parameters;
}

bool GaffOpenMMXMLFileParser::getDihedraLikeData(std::vector<double>& forceConstants, std::vector<double>& phases,
                                                 std::vector<unsigned int>& periodicity, const pugi::xml_node& forceNode) {
  for (pugi::xml_attribute attr : forceNode.attributes()) {
    const std::string name = attr.name();
    if (name.find("k") != std::string::npos && name.size() == 2) {
      forceConstants.push_back(std::stod(attr.value()) * kcalPerMol_per_kJPerMol);
    }
    if (name.find("phase") != std::string::npos) {
      phases.push_back(std::stod(attr.value()) * Utils::Constants::degree_per_rad);
    }
    if (name.find("periodicity") != std::string::npos) {
      periodicity.push_back(std::stoi(attr.value()));
    }
  }
  return forceConstants.size() > 0 && forceConstants.size() == phases.size() && forceConstants.size() == periodicity.size();
}

GaffParameters GaffOpenMMXMLFileParser::parsePeriodicTorsionForce(const pugi::xml_document& xmlDocument, Core::Log& log) {
  GaffParameters parameters;
  pugi::xml_node forces = xmlDocument.child("ForceField").child("PeriodicTorsionForce");
  for (pugi::xml_node force = forces.child("Proper"); force; force = force.next_sibling("Proper")) {
    std::string type1, class1;
    std::string type2, class2;
    std::string type3, class3;
    std::string type4, class4;
    std::vector<double> forceConstants;
    std::vector<double> phases;
    std::vector<unsigned int> periodicities;
    bool success = false;
    try {
      type1 = force.attribute("type1").as_string();
      type2 = force.attribute("type2").as_string();
      type3 = force.attribute("type3").as_string();
      type4 = force.attribute("type4").as_string();
      class1 = force.attribute("class1").as_string();
      class2 = force.attribute("class2").as_string();
      class3 = force.attribute("class3").as_string();
      class4 = force.attribute("class4").as_string();
      success = getDihedraLikeData(forceConstants, phases, periodicities, force);
    }
    catch (...) {
      throw std::runtime_error("Failed to read periodic torsion force definition in XML file.");
    }
    if (!success) {
      throw std::runtime_error("Could not interpret force constants in periodic torsion force data. ");
    }
    if (!class1.empty() || !class2.empty() || !class3.empty() || !class4.empty()) {
      addDihedraTerm(class1, class2, class3, class4, parameters, forceConstants, phases, periodicities, false);
    }
    else {
      if (type1.empty() && type2.empty() && type3.empty() && type4.empty()) {
        throw std::runtime_error("Periodic torsion force information: Missing atom types!");
      }
      addDihedraTerm(type1, type2, type3, type4, parameters, forceConstants, phases, periodicities, false);
    }
  }
  if (parameters.empty()) {
    log.warning << "Warning: "
                << "No (proper) PeriodicTorsionForces were read from the XML file!" << Core::Log::endl;
  }
  return parameters;
}

GaffParameters GaffOpenMMXMLFileParser::parseImproperPeriodicTorsionForce(const pugi::xml_document& xmlDocument,
                                                                          Core::Log& log) {
  GaffParameters parameters;
  pugi::xml_node forces = xmlDocument.child("ForceField").child("PeriodicTorsionForce");
  for (pugi::xml_node force = forces.child("Improper"); force; force = force.next_sibling("Improper")) {
    std::string type1, class1;
    std::string type2, class2;
    std::string type3, class3;
    std::string type4, class4;
    std::vector<double> forceConstants;
    std::vector<double> phases;
    std::vector<unsigned int> periodicities;
    bool success = false;
    try {
      type1 = force.attribute("type1").as_string();
      type2 = force.attribute("type2").as_string();
      type3 = force.attribute("type3").as_string();
      type4 = force.attribute("type4").as_string();
      class1 = force.attribute("class1").as_string();
      class2 = force.attribute("class2").as_string();
      class3 = force.attribute("class3").as_string();
      class4 = force.attribute("class4").as_string();
      success = getDihedraLikeData(forceConstants, phases, periodicities, force);
    }
    catch (...) {
      throw std::runtime_error("Failed read improper periodic torsion bond force definition in XML file.");
    }
    if (!success) {
      throw std::runtime_error("Could not interpret force constants in improper periodic torsion force data. ");
    }
    if (!class1.empty() || !class2.empty() || !class3.empty() || !class4.empty()) {
      addDihedraTerm(class1, class2, class3, class4, parameters, forceConstants, phases, periodicities, true);
    }
    else {
      if (type1.empty() && type2.empty() && type3.empty() && type4.empty()) {
        throw std::runtime_error("Improper periodic torsion force information: Missing atom types!");
      }
      addDihedraTerm(type1, type2, type3, type4, parameters, forceConstants, phases, periodicities, true);
    }
  }
  if (parameters.empty()) {
    log.warning << "Warning: "
                << "No (improper) PeriodicTorsionForces were read from the XML file!" << Core::Log::endl;
  }
  return parameters;
}

GaffParameters GaffOpenMMXMLFileParser::parseNonbondedForce(const pugi::xml_document& xmlDocument, Core::Log& log) {
  GaffParameters parameters;
  pugi::xml_node forces = xmlDocument.child("ForceField").child("NonbondedForce");
  double coulomb14Scale = NAN;
  double lj14Scale = NAN;
  const double sixthRootOfTwoDivTwo = std::pow(2, 1.0 / 6.0) / 2.0;
  try {
    coulomb14Scale = forces.attribute("coulomb14scale").as_double();
    lj14Scale = forces.attribute("lj14scale").as_double();
  }
  catch (...) {
    throw std::runtime_error("Failed to read coulomb13scale or lj14scale definition in XML file.");
  }
  if (std::isnan(coulomb14Scale) || std::isnan(lj14Scale)) {
    throw std::runtime_error("Coulomb or Lennard Jones scaling could not be read from the XML file.");
  }
  parameters.set14ElectrostaticScaling(coulomb14Scale);
  parameters.set14LennardJonesScaling(lj14Scale);
  for (pugi::xml_node force = forces.child("Atom"); force; force = force.next_sibling("Atom")) {
    std::string atomType, atomClass;
    double radius = NAN;
    double epsilon = NAN;
    try {
      atomType = force.attribute("type").as_string();
      atomClass = force.attribute("class").as_string();
      if (force.attribute("sigma")) {
        radius = force.attribute("sigma").as_double() * angstrom_per_nanometer * sixthRootOfTwoDivTwo; // angstrom
      }
      if (force.attribute("epsilon")) {
        epsilon = force.attribute("epsilon").as_double() * kcalPerMol_per_kJPerMol; // kcal/mol
      }
    }
    catch (...) {
      throw std::runtime_error("Failed to read non-bonded force definition in XML file.");
    }
    if (std::isnan(radius) || std::isnan(epsilon)) {
      throw std::runtime_error(
          std::string("Nonbonded force information in XML file incomplete. The following was read from the file:\n") +
          " type=" + atomType + " class " + atomClass + " | radius=" + std::to_string(radius) +
          " | epsilon=" + std::to_string(epsilon));
    }
    if (!atomClass.empty()) {
      parameters.addLennardJones(atomClass, LennardJonesParameters(radius, epsilon));
    }
    else {
      if (atomType.empty()) {
        throw std::runtime_error("Nonbonded force information: Missing atom type!");
      }
      parameters.addLennardJones(atomType, LennardJonesParameters(radius, epsilon));
    }
  }
  if (parameters.empty()) {
    log.warning << "Warning: "
                << "No NonbondedForces were read from the XML file!" << Core::Log::endl;
  }
  return parameters;
}

std::map<std::string, std::string> GaffOpenMMXMLFileParser::parseType2AtomClassMap(const pugi::xml_document& xmlDocument) {
  pugi::xml_node atomTypes = xmlDocument.child("ForceField").child("AtomTypes");
  std::map<std::string, std::string> type2Class;
  for (pugi::xml_node atomType = atomTypes.child("Type"); atomType; atomType = atomType.next_sibling("Type")) {
    const std::string name = atomType.attribute("name").as_string();
    const std::string className = atomType.attribute("class").as_string();
    if (name.empty() || className.empty()) {
      continue;
    }
    type2Class.emplace(name, className);
  }
  return type2Class;
}
void GaffOpenMMXMLFileParser::addDihedraTerm(std::string type1, std::string type2, std::string type3, std::string type4,
                                             GaffParameters& parameters, std::vector<double>& forceConstants,
                                             std::vector<double>& phases, std::vector<unsigned int>& periodicities,
                                             bool improper) {
  type1 = (type1.empty()) ? "X" : type1;
  type2 = (type2.empty()) ? "X" : type2;
  type3 = (type3.empty()) ? "X" : type3;
  type4 = (type4.empty()) ? "X" : type4;
  // Note that it is already ensured that the vectors forceConstants, phases, and periodicities are of the same size.
  for (unsigned int i = 0; i < forceConstants.size(); ++i) {
    if (improper) {
      parameters.addImproperDihedral(ImproperDihedralType(type1, type2, type3, type4),
                                     DihedralParameters(forceConstants[i], phases[i], int(periodicities[i])));
    }
    else {
      parameters.addDihedral(DihedralType(type1, type2, type3, type4),
                             DihedralParameters(forceConstants[i], phases[i], int(periodicities[i])));
    }
  }
}

} // namespace MolecularMechanics
} // namespace Scine
