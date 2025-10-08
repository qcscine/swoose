/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h" // TODO: use this!
#include "Swoose/MolecularMechanics/GAFF/GaffOpenMMXMLFileParser.h"
#include "Swoose/MolecularMechanics/GAFF/GaffParameters.h"
#include <Core/Log.h>
#include <Utils/Constants.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace MolecularMechanics {
namespace Tests {
/**
 * @brief Test the parsing of OpenMM's force field files.
 * @test
 */
class GaffOpenMMXmlFileParserTest : public Test {
 public:
  Core::Log silentLogger = Core::Log::silent();

  static constexpr double bohr_per_nanometer = 10.0 * Utils::Constants::bohr_per_angstrom;
  static constexpr double hartreePerBohr2_per_kJPerMolPerNanometer2 =
      Utils::Constants::hartree_per_kJPerMol / bohr_per_nanometer / bohr_per_nanometer;
};

/**
 * @brief Read a valid OpenMM XML file and check if the tabulated parameters are read correctly.
 */
TEST_F(GaffOpenMMXmlFileParserTest, parseValidFile) {
  const double sixthRootOfTwoDivTwo = std::pow(2, 1.0 / 6.0) / 2.0;
  ASSERT_NO_THROW(GaffOpenMMXMLFileParser::parseParameters(std::vector<std::string>({openmm_test_xml_file}), silentLogger));
  GaffParameters parameters = GaffOpenMMXMLFileParser::parseParameters(openmm_test_xml_file, silentLogger);

  ASSERT_TRUE(parameters.getMMBond("protein-C", "protein-C").hasParameters());
  ASSERT_TRUE(parameters.getMMBond("protein-C", "protein-CA").hasParameters());
  ASSERT_NEAR(parameters.getMMBond("protein-C", "protein-C").getEquilibriumDistance(), 0.1525 * bohr_per_nanometer, 1e-5);
  ASSERT_NEAR(parameters.getMMBond("protein-C", "protein-CA").getEquilibriumDistance(), 0.1409 * bohr_per_nanometer, 1e-5);
  ASSERT_NEAR(parameters.getMMBond("protein-C", "protein-C").getForceConstant(),
              259407.99999999994 * hartreePerBohr2_per_kJPerMolPerNanometer2, 1e-4);
  ASSERT_NEAR(parameters.getMMBond("protein-C", "protein-CA").getForceConstant(),
              392459.19999999995 * hartreePerBohr2_per_kJPerMolPerNanometer2, 1e-4);

  ASSERT_TRUE(parameters.getMMAngle("protein-C", "protein-C", "protein-O").hasParameters());
  ASSERT_TRUE(parameters.getMMAngle("protein-C", "protein-C", "protein-OH").hasParameters());
  ASSERT_NEAR(parameters.getMMAngle("protein-C", "protein-C", "protein-O").getForceConstant(),
              669.44 * Utils::Constants::hartree_per_kJPerMol, 1e-5);
  ASSERT_NEAR(parameters.getMMAngle("protein-C", "protein-C", "protein-OH").getForceConstant(),
              669.44 * Utils::Constants::hartree_per_kJPerMol, 1e-5);
  ASSERT_NEAR(parameters.getMMAngle("protein-C", "protein-C", "protein-O").getEquilibriumAngle(), 2.0943951023931953, 1e-6);
  ASSERT_NEAR(parameters.getMMAngle("protein-C", "protein-C", "protein-OH").getEquilibriumAngle(), 2.0943951023931953, 1e-6);

  ASSERT_TRUE(parameters.getMMDihedrals("X", "protein-C", "protein-C", "X")[0].hasParameters());
  ASSERT_TRUE(parameters.getMMDihedrals("X", "protein-C", "protein-CA", "X")[0].hasParameters());
  ASSERT_NEAR(parameters.getMMDihedrals("X", "protein-C", "protein-C", "X")[0].getHalfBarrierHeight(),
              15.167 * Utils::Constants::hartree_per_kJPerMol, 1e-6);
  ASSERT_NEAR(parameters.getMMDihedrals("X", "protein-C", "protein-CA", "X")[0].getHalfBarrierHeight(),
              15.167 * Utils::Constants::hartree_per_kJPerMol, 1e-6);
  ASSERT_NEAR(parameters.getMMDihedrals("X", "protein-C", "protein-C", "X")[0].getPhaseShift(), 3.141592653589793, 1e-6);
  ASSERT_NEAR(parameters.getMMDihedrals("X", "protein-C", "protein-CA", "X")[0].getPhaseShift(), 3.141592653589793, 1e-6);
  ASSERT_EQ(parameters.getMMDihedrals("X", "protein-C", "protein-C", "X")[0].getPeriodicity(), 2);
  ASSERT_EQ(parameters.getMMDihedrals("X", "protein-C", "protein-CA", "X")[0].getPeriodicity(), 2);

  ASSERT_TRUE(parameters.getMMImproperDihedrals("protein-CA", "protein-2C", "protein-CA", "protein-CA")[0].hasParameters());
  ASSERT_NEAR(parameters.getMMImproperDihedrals("protein-CA", "protein-2C", "protein-CA", "protein-CA")[0].getHalfBarrierHeight(),
              4.6024 * Utils::Constants::hartree_per_kJPerMol, 1e-6);
  ASSERT_NEAR(parameters.getMMImproperDihedrals("protein-CA", "protein-2C", "protein-CA", "protein-CA")[0].getPhaseShift(),
              3.141592653589793, 1e-6);
  ASSERT_EQ(parameters.getMMImproperDihedrals("protein-CA", "protein-2C", "protein-CA", "protein-CA")[0].getPeriodicity(), 2);

  ASSERT_NO_THROW(parameters.getMMLennardJones("protein-C", "protein-C"));
  ASSERT_NO_THROW(parameters.getMMLennardJones("protein-C", "protein-C8"));
  ASSERT_NO_THROW(parameters.getMMLennardJones("protein-C8", "protein-C8"));
  const std::unordered_map<std::string, LennardJonesParameters>& ljParams = parameters.getLennardJonesParameters();
  const auto itLjPC = ljParams.find("protein-C");
  const auto itLjPC8 = ljParams.find("protein-C8");
  ASSERT_NEAR(itLjPC->second.getWellDepth(), 0.359824 * Utils::Constants::hartree_per_kJPerMol, 1e-6);
  ASSERT_NEAR(itLjPC->second.getVdwRadius(), 0.3399669508423535 * bohr_per_nanometer * sixthRootOfTwoDivTwo, 1e-6);
  ASSERT_NEAR(itLjPC8->second.getWellDepth(), 0.4577296 * Utils::Constants::hartree_per_kJPerMol, 1e-6);
  ASSERT_NEAR(itLjPC8->second.getVdwRadius(), 0.3399669508423535 * bohr_per_nanometer * sixthRootOfTwoDivTwo, 1e-6);
}

TEST_F(GaffOpenMMXmlFileParserTest, ensureFailureForBrokenFiles) {
  ASSERT_THROW(GaffOpenMMXMLFileParser::parseParameters(openmm_broken_xml_file_bonds, silentLogger), std::runtime_error);
  ASSERT_THROW(GaffOpenMMXMLFileParser::parseParameters(openmm_broken_xml_file_angles, silentLogger), std::runtime_error);
  ASSERT_THROW(GaffOpenMMXMLFileParser::parseParameters(openmm_broken_xml_file_dihedra, silentLogger), std::runtime_error);
  ASSERT_THROW(GaffOpenMMXMLFileParser::parseParameters(openmm_broken_xml_file_improper, silentLogger), std::runtime_error);
  ASSERT_THROW(GaffOpenMMXMLFileParser::parseParameters(openmm_broken_xml_file_lj, silentLogger), std::runtime_error);
}

TEST_F(GaffOpenMMXmlFileParserTest, atomTypeClasses) {
  /*
   * This parameter file is a bit more nasty. Most interactions are encoded by atom classes and not by atom type.
   */
  auto parameters = GaffOpenMMXMLFileParser::parseParameters(openmm_amber99sbildn_xml_file, silentLogger);

  ASSERT_TRUE(parameters.getMMBond("C", "C").hasParameters());
  ASSERT_TRUE(parameters.getMMBond("6", "6").hasParameters());  // Encoded as classes C-C
  ASSERT_TRUE(parameters.getMMBond("6", "33").hasParameters()); // Encoded as classes C-C
  ASSERT_TRUE(parameters.getMMBond("N", "C").hasParameters());
  ASSERT_TRUE(parameters.getMMBond("0", "25").hasParameters()); // Encoded as classes N-C

  ASSERT_TRUE(parameters.getMMAngle("C", "C", "O").hasParameters());
  ASSERT_TRUE(parameters.getMMAngle("6", "25", "7").hasParameters()); // Encoded as classes C-C-O

  ASSERT_NO_THROW(parameters.getMMDihedrals("N", "C6", "CT", "C"));
  ASSERT_NO_THROW(parameters.getMMDihedrals("0", "57", "10", "33")); // Encoded as classes X-C6-CT-X

  ASSERT_NO_THROW(parameters.getMMImproperDihedrals("C5", "N", "C", "O"));
  ASSERT_NO_THROW(parameters.getMMImproperDihedrals("45", "0", "6", "7")); // Encoded as classes C5-X-X-O

  ASSERT_NO_THROW(parameters.getMMLennardJones("0", "1700")); // The Lennard Jones parameters are all encoded by type.
  ASSERT_NO_THROW(parameters.getMMLennardJones("645", "27"));
}

} // namespace Tests
} // namespace MolecularMechanics
} // namespace Scine
