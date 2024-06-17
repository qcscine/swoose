/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include <Core/Log.h>
#include <Swoose/MMParametrization/MMParametrizationSettings.h>
#include <Swoose/MMParametrization/MolecularSystemPartitioner.h>
#include <Swoose/MMParametrization/ParametrizationUtils/ConnectivityGenerator.h>
#include <Swoose/MMParametrization/ParametrizationUtils/ReparametrizationHelper.h>
#include <Swoose/MMParametrization/Parametrizer.h>
#include <Swoose/MolecularMechanics/SFAM/SfamAtomTypeIdentifier.h>
#include <Swoose/MolecularMechanics/SFAM/SfamParameterParser.h>
#include <Swoose/MolecularMechanics/Topology/IndexedStructuralTopologyCreator.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/FilesystemHelpers.h>
#include <Utils/IO/NativeFilenames.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>

using namespace testing;
namespace Scine {
using namespace MMParametrization;
namespace Tests {

/**
 * @class AParametrizationOfSmallSystemsTest ParametrizationOfSmallSystemsTest.cpp
 * @brief Tests the automated parametrization for small molecular systems.
 * @test
 */
class AParametrizationOfSmallSystemsTest : public Test {
 public:
  void SetUp() override {
    boost::filesystem::remove_all(connFile);
  }
  void TearDown() override {
    boost::filesystem::remove_all(connFile);
    boost::filesystem::remove_all(parFile);
  }
  std::string connFile = "Connectivity.dat";
  std::string parFile = "Parameters.dat";
  Parametrizer parametrizer;
  Core::Log silentLogger = Core::Log::silent();
};

TEST_F(AParametrizationOfSmallSystemsTest, ReferenceDataFilesAreCorrectlyGeneratedForAlanine) {
  parametrizer.setLog(silentLogger);
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;

  parametrizer.settings().modifyInt(SwooseUtilities::SettingsNames::numberAtomsThreshold, 120);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, alanine_ref_calc_dir);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "write");
  parametrizer.parametrize(structure);

  // Check that the info file was correctly written.
  std::string line;
  std::ifstream infoFile;
  infoFile.open(Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "0", "info.dat"));
  if (infoFile.is_open()) {
    getline(infoFile, line);
  }
  infoFile.close();
  ASSERT_THAT(line, Eq("0  1"));

  // Check that the molecular structure was correctly written.
  auto s = Utils::ChemicalFileHandler::read(Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "0", "molecule.xyz"))
               .first;

  ASSERT_THAT(s.size(), Eq(structure.size()));

  Utils::PositionCollection p1 = s.getPositions();
  Utils::PositionCollection p2 = structure.getPositions();
  for (int i = 0; i < p1.rows(); ++i) {
    ASSERT_TRUE(Utils::ElementInfo::Z(s.getElement(i)) == Utils::ElementInfo::Z(structure.getElement(i)));
    for (int j = 0; j < p1.cols(); ++j) {
      ASSERT_THAT(p1(i, j), DoubleNear(p2(i, j), 1e-4));
    }
  }

  // Delete the generated files
  boost::filesystem::remove_all(Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "0", "info.dat"));
  boost::filesystem::remove_all(Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "0", "molecule.xyz"));
}

TEST_F(AParametrizationOfSmallSystemsTest, AdditionalReferenceDataAreWrittenForTitrableSites) {
  parametrizer.setLog(silentLogger);
  // Glutamine is an acid, so the non-ref form is deprotonated
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(glutamine_xyz_file).first;

  parametrizer.settings().modifyInt(SwooseUtilities::SettingsNames::numberAtomsThreshold, 120);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, glutamine_ref_calc_dir);
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::titrate, true);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "write");

  // first, parametrizer should throw an error when no titration site file is provided
  ASSERT_THROW(parametrizer.parametrize(structure), std::runtime_error);

  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::titrationSiteFile, glutamine_titration_site_file);
  parametrizer.parametrize(structure);
  // check that an additional directory was created for the non-ref state
  auto refDir = Utils::NativeFilenames::combinePathSegments(glutamine_ref_calc_dir, "0");
  auto nonRefDir = Utils::NativeFilenames::combinePathSegments(glutamine_ref_calc_dir, "0", "non_ref_state");
  ASSERT_TRUE(boost::filesystem::exists(refDir));
  ASSERT_TRUE(boost::filesystem::exists(nonRefDir));

  ASSERT_TRUE(boost::filesystem::exists(Utils::NativeFilenames::combinePathSegments(refDir, "molecule.xyz")));
  ASSERT_TRUE(boost::filesystem::exists(Utils::NativeFilenames::combinePathSegments(refDir, "info.dat")));

  ASSERT_TRUE(boost::filesystem::exists(Utils::NativeFilenames::combinePathSegments(nonRefDir, "molecule.xyz")));
  ASSERT_TRUE(boost::filesystem::exists(Utils::NativeFilenames::combinePathSegments(nonRefDir, "info.dat")));

  // Check that the info file was correctly written.
  std::string line;
  std::ifstream infoFile;
  infoFile.open(Utils::NativeFilenames::combinePathSegments(refDir, "info.dat"));
  if (infoFile.is_open()) {
    getline(infoFile, line);
  }
  infoFile.close();
  ASSERT_THAT(line, Eq("0  1"));

  // Check that the info file was correctly written for the non-ref structure.
  std::ifstream infoFileNonRef;
  infoFileNonRef.open(Utils::NativeFilenames::combinePathSegments(nonRefDir, "info.dat"));
  if (infoFileNonRef.is_open()) {
    getline(infoFileNonRef, line);
  }
  infoFileNonRef.close();
  ASSERT_THAT(line, Eq("-1  1"));

  // Check that the molecular structure was correctly written for the non-ref state
  auto deprotStruct =
      Utils::ChemicalFileHandler::read(Utils::NativeFilenames::combinePathSegments(nonRefDir, "molecule.xyz")).first;

  // Check that the molecular structure was correctly written for the non-ref state
  auto protStruct =
      Utils::ChemicalFileHandler::read(Utils::NativeFilenames::combinePathSegments(refDir, "molecule.xyz")).first;

  ASSERT_THAT(deprotStruct.size(), Eq(protStruct.size() - 1));

  // Delete the generated files
  boost::filesystem::remove_all(Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "0", "info.dat"));
  boost::filesystem::remove_all(Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "0", "molecule.xyz"));
}

TEST_F(AParametrizationOfSmallSystemsTest, GlutamineIsCorrectlyParametrizedIncludingTitration) {
  parametrizer.setLog(silentLogger);
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(glutamine_xyz_file).first;

  parametrizer.settings().modifyInt(SwooseUtilities::SettingsNames::numberAtomsThreshold, 120);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, glutamine_ref_calc_dir);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  parametrizer.settings().modifyString(Utils::SettingsNames::parameterFilePath, parFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "direct");
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useGaussianOptionKey, false);
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useCsvInputFormat, false);
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::titrate, true);
  // parametrization including titration does not work in direct mode
  ASSERT_THROW(parametrizer.parametrize(structure), std::runtime_error);

  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  // TODO
}

TEST_F(AParametrizationOfSmallSystemsTest, AlanineIsCorrectlyParametrizedWithGaussian) {
  parametrizer.setLog(silentLogger);

  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;

  parametrizer.settings().modifyInt(SwooseUtilities::SettingsNames::numberAtomsThreshold, 120);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, alanine_ref_calc_dir);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  parametrizer.settings().modifyString(Utils::SettingsNames::parameterFilePath, parFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useGaussianOptionKey, true);
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useCsvInputFormat, false);
  parametrizer.parametrize(structure);

  // Check connectivity file
  auto listsOfNeighbors = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(connFile);
  for (const auto& atom : listsOfNeighbors) {
    ASSERT_FALSE(atom.empty());
  }

  // Check that the parameter file was correctly written.
  std::string line;
  std::ifstream paramFile;
  paramFile.open(parFile);
  if (paramFile.is_open()) {
    getline(paramFile, line);
  }
  paramFile.close();
  std::string infoString = "# MM parameters of SFAM generated by SCINE";
  ASSERT_THAT(line, Eq(infoString));

  // Generate topology and atom types
  MolecularMechanics::IndexedStructuralTopologyCreator topologyCreator(listsOfNeighbors);
  auto topology = topologyCreator.calculateIndexedStructuralTopology();
  topologyCreator.addHydrogenBondsToIndexedStructuralTopology(topology, structure);
  MolecularMechanics::SfamAtomTypeIdentifier atomTypeIdentifier(structure.size(), structure.getElements(), listsOfNeighbors);
  MolecularMechanics::SfamAtomTypeLevel atl =
      MolecularMechanics::SfamAtomTypeIdentifier::generateSfamAtomTypeLevelFromString("high");
  auto atomTypes = atomTypeIdentifier.getAtomTypes(atl);
  ASSERT_THAT(atomTypes.size(), Eq(structure.size()));

  MolecularMechanics::SfamParameterParser parser(parFile, atomTypes);
  auto parameters = parser.parseParameters();

  ASSERT_THAT(parameters->getNonCovalentParameters().size(), Eq(5));

  // Check that the bond parameters are existing and have reasonable values
  for (const auto& bond : parameters->getBonds()) {
    bool reasonableForceConstant = (bond.second.getForceConstant() > 200) && (bond.second.getForceConstant() < 2000);
    ASSERT_TRUE(reasonableForceConstant);
    bool reasonableEqBondDistance =
        (bond.second.getEquilibriumBondLength() > 0.5) && (bond.second.getEquilibriumBondLength() < 3.0);
    ASSERT_TRUE(reasonableEqBondDistance);
  }

  // Check that the angle parameters are existing and have reasonable values
  for (const auto& angle : parameters->getAngles()) {
    bool reasonableForceConstant = (angle.second.getForceConstant() > 70) && (angle.second.getForceConstant() < 500);
    ASSERT_TRUE(reasonableForceConstant);
    bool reasonableEqAngle = (angle.second.getEquilibriumAngle() >= 0) && (angle.second.getEquilibriumAngle() <= 180);
    ASSERT_TRUE(reasonableEqAngle);
  }

  // Check that the dihedral parameters are existing and have reasonable values
  for (const auto& dihedral : parameters->getDihedrals()) {
    bool reasonableHalfBarrierHeight = std::abs(dihedral.second.getHalfBarrierHeight()) < 7;
    ASSERT_TRUE(reasonableHalfBarrierHeight);
    bool reasonablePhaseShift = (dihedral.second.getPhaseShift() >= 0) && (dihedral.second.getPhaseShift() <= 180);
    ASSERT_TRUE(reasonablePhaseShift);
    bool reasonablePeriodicity = (dihedral.second.getPeriodicity() > 0) && (dihedral.second.getPeriodicity() < 10);
    ASSERT_TRUE(reasonablePeriodicity);
  }

  ASSERT_FALSE(parameters->getImproperDihedrals().empty());

  // Check atomic charges and c6 coefficients
  double sumOfAtomicCharges = 0.0;
  std::vector<double> atomicCharges = parameters->getChargesForEachAtom(atomTypes);
  ASSERT_THAT(atomicCharges.size(), Eq(structure.size()));
  for (int i = 0; i < structure.size(); ++i) {
    double charge = atomicCharges[i];
    ASSERT_TRUE(std::abs(charge) < 1);
    sumOfAtomicCharges += charge;
    for (int j = 0; j < i; ++j) {
      ASSERT_TRUE(parameters->getC6(atomTypes.getAtomType(i), atomTypes.getAtomType(j)) > 0);
    }
  }

  // Sum of atomic charges should be roughly zero
  ASSERT_TRUE(std::abs(sumOfAtomicCharges) < 1e-5);

  // Copy needed for test below
  std::string existingParFile = "existing_parameters.dat";
  Utils::FilesystemHelpers::copyFile(parFile, existingParFile);

  // Delete the generated files
  boost::filesystem::remove_all(connFile);
  boost::filesystem::remove_all(parFile);

  // Redo parametrization, but with re-parametrization mode:
  Parametrizer reparametrizer;
  reparametrizer.setLog(silentLogger);
  reparametrizer.settings().modifyString(SwooseUtilities::SettingsNames::existingParameters, existingParFile);
  reparametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, alanine_ref_calc_dir);
  reparametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  reparametrizer.settings().modifyString(Utils::SettingsNames::parameterFilePath, parFile);
  reparametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  reparametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useGaussianOptionKey, false);
  reparametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useCsvInputFormat, false);
  reparametrizer.parametrize(structure);

  ASSERT_TRUE(boost::filesystem::exists(connFile));
  ASSERT_TRUE(boost::filesystem::exists(parFile));
  ASSERT_TRUE(boost::filesystem::exists(existingParFile));

  // Delete the generated existing parameters file
  boost::filesystem::remove_all(existingParFile);
}

TEST_F(AParametrizationOfSmallSystemsTest, AlanineIsCorrectlyParametrizedWithoutGaussian) {
  parametrizer.setLog(silentLogger);
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;

  parametrizer.settings().modifyInt(SwooseUtilities::SettingsNames::numberAtomsThreshold, 120);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, alanine_ref_calc_dir);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  parametrizer.settings().modifyString(Utils::SettingsNames::parameterFilePath, parFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useGaussianOptionKey, false);
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useCsvInputFormat, false);
  parametrizer.parametrize(structure);

  // Check connectivity file
  auto listsOfNeighbors = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(connFile);
  for (const auto& atom : listsOfNeighbors) {
    ASSERT_FALSE(atom.empty());
  }

  // Generate topology and atom types
  MolecularMechanics::IndexedStructuralTopologyCreator topologyCreator(listsOfNeighbors);
  auto topology = topologyCreator.calculateIndexedStructuralTopology();
  topologyCreator.addHydrogenBondsToIndexedStructuralTopology(topology, structure);
  MolecularMechanics::SfamAtomTypeIdentifier atomTypeIdentifier(structure.size(), structure.getElements(), listsOfNeighbors);
  MolecularMechanics::SfamAtomTypeLevel atl =
      MolecularMechanics::SfamAtomTypeIdentifier::generateSfamAtomTypeLevelFromString("high");
  auto atomTypes = atomTypeIdentifier.getAtomTypes(atl);
  ASSERT_THAT(atomTypes.size(), Eq(structure.size()));

  MolecularMechanics::SfamParameterParser parser(parFile, atomTypes);
  auto parameters = parser.parseParameters();

  ASSERT_THAT(parameters->getNonCovalentParameters().size(), Eq(5));

  // Check that the bond parameters are existing and have reasonable values
  for (const auto& bond : parameters->getBonds()) {
    bool reasonableForceConstant = (bond.second.getForceConstant() > 200) && (bond.second.getForceConstant() < 2000);
    ASSERT_TRUE(reasonableForceConstant);
    bool reasonableEqBondDistance =
        (bond.second.getEquilibriumBondLength() > 0.5) && (bond.second.getEquilibriumBondLength() < 3.0);
    ASSERT_TRUE(reasonableEqBondDistance);
  }

  // Check that the angle parameters are existing and have reasonable values
  for (const auto& angle : parameters->getAngles()) {
    bool reasonableForceConstant = (angle.second.getForceConstant() > 70) && (angle.second.getForceConstant() < 500);
    ASSERT_TRUE(reasonableForceConstant);
    bool reasonableEqAngle = (angle.second.getEquilibriumAngle() >= 0) && (angle.second.getEquilibriumAngle() <= 180);
    ASSERT_TRUE(reasonableEqAngle);
  }

  // Check that the dihedral parameters are existing and have reasonable values
  for (const auto& dihedral : parameters->getDihedrals()) {
    bool reasonableHalfBarrierHeight = std::abs(dihedral.second.getHalfBarrierHeight()) < 7;
    ASSERT_TRUE(reasonableHalfBarrierHeight);
    bool reasonablePhaseShift = (dihedral.second.getPhaseShift() >= 0) && (dihedral.second.getPhaseShift() <= 180);
    ASSERT_TRUE(reasonablePhaseShift);
    bool reasonablePeriodicity = (dihedral.second.getPeriodicity() > 0) && (dihedral.second.getPeriodicity() < 10);
    ASSERT_TRUE(reasonablePeriodicity);
  }

  ASSERT_FALSE(parameters->getImproperDihedrals().empty());

  // Check atomic charges and c6 coefficients
  double sumOfAtomicCharges = 0.0;
  std::vector<double> atomicCharges = parameters->getChargesForEachAtom(atomTypes);
  for (int i = 0; i < structure.size(); ++i) {
    double charge = atomicCharges[i];
    ASSERT_TRUE(std::abs(charge) < 1);
    sumOfAtomicCharges += charge;
    for (int j = 0; j < i; ++j) {
      ASSERT_TRUE(parameters->getC6(atomTypes.getAtomType(i), atomTypes.getAtomType(j)) > 0);
    }
  }

  // Sum of atomic charges should be roughly zero
  ASSERT_TRUE(std::abs(sumOfAtomicCharges) < 1e-5);

  // Copy needed for test below
  std::string existingParFile = "existing_parameters.dat";
  Utils::FilesystemHelpers::copyFile(parFile, existingParFile);

  // Delete the generated files
  boost::filesystem::remove_all(connFile);
  boost::filesystem::remove_all(parFile);

  // Redo parametrization, but with re-parametrization mode:
  Parametrizer reparametrizer;
  reparametrizer.setLog(silentLogger);
  reparametrizer.settings().modifyString(SwooseUtilities::SettingsNames::existingParameters, existingParFile);
  reparametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, alanine_ref_calc_dir);
  reparametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  reparametrizer.settings().modifyString(Utils::SettingsNames::parameterFilePath, parFile);
  reparametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  reparametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useGaussianOptionKey, false);
  reparametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useCsvInputFormat, false);
  reparametrizer.parametrize(structure);

  ASSERT_TRUE(boost::filesystem::exists(connFile));
  ASSERT_TRUE(boost::filesystem::exists(parFile));
  ASSERT_TRUE(boost::filesystem::exists(existingParFile));

  // Delete the generated existing parameters file
  boost::filesystem::remove_all(existingParFile);
}

TEST_F(AParametrizationOfSmallSystemsTest, AlanineIsCorrectlyParametrizedFromDataInCsvFiles) {
  parametrizer.setLog(silentLogger);
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;

  parametrizer.settings().modifyInt(SwooseUtilities::SettingsNames::numberAtomsThreshold, 120);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, alanine_ref_calc_dir);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  parametrizer.settings().modifyString(Utils::SettingsNames::parameterFilePath, parFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useGaussianOptionKey, false);
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useCsvInputFormat, true);
  parametrizer.parametrize(structure);

  // Check connectivity file
  auto listsOfNeighbors = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(connFile);
  for (const auto& atom : listsOfNeighbors) {
    ASSERT_FALSE(atom.empty());
  }

  // Generate topology and atom types
  MolecularMechanics::IndexedStructuralTopologyCreator topologyCreator(listsOfNeighbors);
  auto topology = topologyCreator.calculateIndexedStructuralTopology();
  topologyCreator.addHydrogenBondsToIndexedStructuralTopology(topology, structure);
  MolecularMechanics::SfamAtomTypeIdentifier atomTypeIdentifier(structure.size(), structure.getElements(), listsOfNeighbors);
  MolecularMechanics::SfamAtomTypeLevel atl =
      MolecularMechanics::SfamAtomTypeIdentifier::generateSfamAtomTypeLevelFromString("high");
  auto atomTypes = atomTypeIdentifier.getAtomTypes(atl);
  ASSERT_THAT(atomTypes.size(), Eq(structure.size()));

  MolecularMechanics::SfamParameterParser parser(parFile, atomTypes);
  auto parameters = parser.parseParameters();

  ASSERT_THAT(parameters->getNonCovalentParameters().size(), Eq(5));

  // Check that the bond parameters are existing and have reasonable values
  for (const auto& bond : parameters->getBonds()) {
    bool reasonableForceConstant = (bond.second.getForceConstant() > 200) && (bond.second.getForceConstant() < 2000);
    ASSERT_TRUE(reasonableForceConstant);
    bool reasonableEqBondDistance =
        (bond.second.getEquilibriumBondLength() > 0.5) && (bond.second.getEquilibriumBondLength() < 3.0);
    ASSERT_TRUE(reasonableEqBondDistance);
  }

  // Check that the angle parameters are existing and have reasonable values
  for (const auto& angle : parameters->getAngles()) {
    bool reasonableForceConstant = (angle.second.getForceConstant() > 70) && (angle.second.getForceConstant() < 500);
    ASSERT_TRUE(reasonableForceConstant);
    bool reasonableEqAngle = (angle.second.getEquilibriumAngle() >= 0) && (angle.second.getEquilibriumAngle() <= 180);
    ASSERT_TRUE(reasonableEqAngle);
  }

  // Check that the dihedral parameters are existing and have reasonable values
  for (const auto& dihedral : parameters->getDihedrals()) {
    bool reasonableHalfBarrierHeight = std::abs(dihedral.second.getHalfBarrierHeight()) < 7;
    ASSERT_TRUE(reasonableHalfBarrierHeight);
    bool reasonablePhaseShift = (dihedral.second.getPhaseShift() >= 0) && (dihedral.second.getPhaseShift() <= 180);
    ASSERT_TRUE(reasonablePhaseShift);
    bool reasonablePeriodicity = (dihedral.second.getPeriodicity() > 0) && (dihedral.second.getPeriodicity() < 10);
    ASSERT_TRUE(reasonablePeriodicity);
  }

  ASSERT_FALSE(parameters->getImproperDihedrals().empty());

  // Check atomic charges and c6 coefficients
  double sumOfAtomicCharges = 0.0;
  std::vector<double> atomicCharges = parameters->getChargesForEachAtom(atomTypes);
  for (int i = 0; i < structure.size(); ++i) {
    double charge = atomicCharges[i];
    ASSERT_TRUE(std::abs(charge) < 1);
    sumOfAtomicCharges += charge;
    for (int j = 0; j < i; ++j) {
      ASSERT_TRUE(parameters->getC6(atomTypes.getAtomType(i), atomTypes.getAtomType(j)) > 0);
    }
  }

  // Sum of atomic charges should be roughly zero
  ASSERT_TRUE(std::abs(sumOfAtomicCharges) < 1e-5);

  // Copy needed for test below
  std::string existingParFile = "existing_parameters.dat";
  Utils::FilesystemHelpers::copyFile(parFile, existingParFile);

  // Delete the generated files
  boost::filesystem::remove_all(connFile);
  boost::filesystem::remove_all(parFile);
}

// Run this test only in release builds
#ifdef NDEBUG
TEST_F(AParametrizationOfSmallSystemsTest, EquilibriumValuesAreCorrectlyAssignedForImproperDihedrals) {
  parametrizer.setLog(silentLogger);
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(guanidine_fragment_xyz_file).first;

  auto infoFile = Utils::NativeFilenames::combinePathSegments(guanidine_fragment_ref_calc_dir, "0", "atomic_info.dat");
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, guanidine_fragment_ref_calc_dir);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  parametrizer.settings().modifyString(Utils::SettingsNames::parameterFilePath, parFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::atomicInformationFile, infoFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useGaussianOptionKey, false);
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useCsvInputFormat, false);
  parametrizer.parametrize(structure);

  // Check connectivity file
  auto listsOfNeighbors = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(connFile);
  for (const auto& atom : listsOfNeighbors) {
    ASSERT_FALSE(atom.empty());
  }

  // Generate topology and atom types
  MolecularMechanics::IndexedStructuralTopologyCreator topologyCreator(listsOfNeighbors);
  auto topology = topologyCreator.calculateIndexedStructuralTopology();
  topologyCreator.addHydrogenBondsToIndexedStructuralTopology(topology, structure);
  MolecularMechanics::SfamAtomTypeIdentifier atomTypeIdentifier(structure.size(), structure.getElements(), listsOfNeighbors);
  MolecularMechanics::SfamAtomTypeLevel atl =
      MolecularMechanics::SfamAtomTypeIdentifier::generateSfamAtomTypeLevelFromString("high");
  auto atomTypes = atomTypeIdentifier.getAtomTypes(atl);
  ASSERT_THAT(atomTypes.size(), Eq(structure.size()));

  MolecularMechanics::SfamParameterParser parser(parFile, atomTypes);
  auto parameters = parser.parseParameters();

  // Check that the improper dihedral parameters are existing and have reasonable values
  for (const auto& improper_dihedral : parameters->getImproperDihedrals())
    ASSERT_TRUE(improper_dihedral.second.getEquilibriumAngle() == 0);
}
#endif

TEST_F(AParametrizationOfSmallSystemsTest, CustomConnectivityIsImportedCorrectly) {
  // parametrizer.setLog(silentLogger);
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(guanidine_fragment_xyz_file).first;

  // Write custom connectivity file
  std::ofstream customConnectivityFile(connFile);
  for (int i = 0; i < structure.size(); ++i) {
    if (i == 0) {
      customConnectivityFile << "1\n";
      continue;
    }
    if (i == 1) {
      customConnectivityFile << "0\n";
      continue;
    }
    customConnectivityFile << "-1\n";
  }
  customConnectivityFile.close();

  auto infoFile = Utils::NativeFilenames::combinePathSegments(guanidine_fragment_ref_calc_dir, "0", "atomic_info.dat");
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, guanidine_fragment_ref_calc_dir);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  parametrizer.settings().modifyString(Utils::SettingsNames::parameterFilePath, parFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::atomicInformationFile, infoFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useGaussianOptionKey, false);
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::refineConnectivity, false); // Important for this
                                                                                                 // test
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useCsvInputFormat, false);
  parametrizer.parametrize(structure);

  // Read connectivity from file
  auto listsOfNeighbors = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(connFile);

  // Generate topology and atom types
  MolecularMechanics::IndexedStructuralTopologyCreator topologyCreator(listsOfNeighbors);
  auto topology = topologyCreator.calculateIndexedStructuralTopology();
  topologyCreator.addHydrogenBondsToIndexedStructuralTopology(topology, structure);
  MolecularMechanics::SfamAtomTypeIdentifier atomTypeIdentifier(structure.size(), structure.getElements(), listsOfNeighbors);
  MolecularMechanics::SfamAtomTypeLevel atl =
      MolecularMechanics::SfamAtomTypeIdentifier::generateSfamAtomTypeLevelFromString("high");
  auto atomTypes = atomTypeIdentifier.getAtomTypes(atl);
  ASSERT_THAT(atomTypes.size(), Eq(structure.size()));

  MolecularMechanics::SfamParameterParser parser(parFile, atomTypes);
  auto parameters = parser.parseParameters();

  // Check that the parameters contain only one bond
  ASSERT_THAT(parameters->getBonds().size(), Eq(1)); // only one bond was specified in the custom conn. file above
}

TEST_F(AParametrizationOfSmallSystemsTest, ConnectivityFileWithSelfBondedAtomsCausesException) {
  parametrizer.setLog(silentLogger);
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(guanidine_fragment_xyz_file).first;

  // Write custom connectivity file
  std::ofstream customConnectivityFile(connFile);
  for (int i = 0; i < structure.size(); ++i) {
    if (i == 0) {
      customConnectivityFile << "1\n";
      continue;
    }
    if (i == 1) {
      customConnectivityFile << "0\n";
      continue;
    }
    if (i == 10) {
      customConnectivityFile << "10 7 5\n";
      continue;
    }
    customConnectivityFile << "-1\n";
  }
  customConnectivityFile.close();

  auto infoFile = Utils::NativeFilenames::combinePathSegments(guanidine_fragment_ref_calc_dir, "0", "atomic_info.dat");
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataDirectory, guanidine_fragment_ref_calc_dir);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);
  parametrizer.settings().modifyString(Utils::SettingsNames::parameterFilePath, parFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::atomicInformationFile, infoFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode, "read");
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::useGaussianOptionKey, false);

  std::string exceptionString = "";
  try {
    parametrizer.parametrize(structure);
  }
  catch (const std::runtime_error& e) {
    exceptionString = e.what();
  }

  ASSERT_STREQ(exceptionString.c_str(), "Error in connectivity file. Atom 10 is bonded to itself.");
}

} // namespace Tests
} // namespace Scine
