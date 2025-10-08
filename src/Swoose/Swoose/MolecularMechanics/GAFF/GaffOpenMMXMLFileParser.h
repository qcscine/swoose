/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_GAFFOPENMMXMLFILEPARSER_H
#define SWOOSE_GAFFOPENMMXMLFILEPARSER_H

#include <Core/Log.h>
#include <map>
#include <pugixml.hpp>
#include <string>
#include <vector>

namespace Scine {

namespace MolecularMechanics {

class GaffParameters;

/**
 * @class GaffOpenMMXMLFileParser GaffOpenMMXMLFileParser.h
 * @brief Read a force field from an OpenMM XML file. At the moment this is tailored to Amber force fields.
 *
 * The following OpenMM force types are supported:
 * HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce, NonbondedForce.
 *
 * A valid XML file may look like this:
 * <ForceField>
 *   <HarmonicBondForce>
 *     <Bond k="259407.99999999994" length="0.1525" type1="protein-C" type2="protein-C"/>
 *   </HarmonicBondForce>
 *   <NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5">
 *     <Atom epsilon="0.0" sigma="1.0" type="protein-HO"/>
 *   </NonbondedForce>
 * </ForceField>
 *
 * Any residue or atom type definitions in the file are ignored.
 */
class GaffOpenMMXMLFileParser {
 private:
  // Purely static.
  GaffOpenMMXMLFileParser() = default;
  ~GaffOpenMMXMLFileParser() = default;

 public:
  /**
   * @brief Read the parameters from a single XML file.
   * @param xmlFilePath The XML file path.
   * @param log         The log to issue warnings.
   * @return The parameters.
   */
  static GaffParameters parseParameters(const std::string& xmlFilePath, Core::Log& log);
  /**
   * @brief Read the parameters from multiple XML files. The parameters are not checked for duplicated force
   *   definitions.
   * @param xmlFilePaths The XML file paths.
   * @param log          The log to issue warnings.
   * @return The parameters.
   */
  static GaffParameters parseParameters(const std::vector<std::string>& xmlFilePaths, Core::Log& log);

 private:
  // Each of these functions reads a distinct type of forces, e.g., harmonic bond forces.
  static GaffParameters parseHarmonicBondForce(const pugi::xml_document& xmlDocument, Core::Log& log);
  static GaffParameters parseHarmonicAngleForce(const pugi::xml_document& xmlDocument, Core::Log& log);
  static GaffParameters parsePeriodicTorsionForce(const pugi::xml_document& xmlDocument, Core::Log& log);
  static GaffParameters parseImproperPeriodicTorsionForce(const pugi::xml_document& xmlDocument, Core::Log& log);
  static GaffParameters parseNonbondedForce(const pugi::xml_document& xmlDocument, Core::Log& log);
  static std::map<std::string, std::string> parseType2AtomClassMap(const pugi::xml_document& xmlDocument);
  static void addDihedraTerm(std::string type1, std::string type2, std::string type3, std::string type4,
                             GaffParameters& parameters, std::vector<double>& forceConstants,
                             std::vector<double>& phases, std::vector<unsigned int>& periodicities, bool improper);
  /**
   * @brief Extract force constants, phases, and periodicities for dihedra angles or improper dihedra angles from an xml
   * force node. Already include the unit conversion from kJ/mol to kcal/mol and rad to degree.
   * @param forceConstants The force constants will be written to this vector.
   * @param phases         The phases will be written to this vector.
   * @param periodicity    The periodicity information will be written to this vector.
   * @param forceNode      The XML node encoding the force information.
   * @return Returns true if force constants were extracted and the information is consistent between phases, force
   * constants, and periodicities. Returns false if something went wrong.
   */
  static bool getDihedraLikeData(std::vector<double>& forceConstants, std::vector<double>& phases,
                                 std::vector<unsigned int>& periodicity, const pugi::xml_node& forceNode);
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // SWOOSE_GAFFOPENMMXMLFILEPARSER_H
