/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include "TestUtilities/AminoAcidDataForTests.h"
#include <Core/Log.h>
#include <Swoose/StructurePreparation/AminoAcidGraphRepresentations.h>
#include <Swoose/StructurePreparation/Protonation/ProtonationHandler.h>
#include <Swoose/StructurePreparation/Protonation/TitrationData.h>
#include <Swoose/StructurePreparation/StructurePreparationData.h>
#include <Swoose/StructurePreparation/StructurePreparationHelper.h>
#include <Swoose/StructurePreparation/StructurePreparationIO.h>
#include <Swoose/StructurePreparation/StructurePreparationSettings.h>
#include <Swoose/StructurePreparation/StructureProcessor.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/SettingsNames.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/ChemicalFileFormats/OpenBabelStreamHandler.h>
#include <Utils/IO/ChemicalFileFormats/PdbStreamHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>
#include <string>

using namespace testing;
namespace Scine {
using namespace StructurePreparation;
namespace Tests {

/**
 * @class AStructurePreparationTest PdbPreparationTest.cpp
 * @class Tests the automated PDB preparation for plastocyanin.
 */

class AStructurePreparationTest : public Test {
 public:
  StructureProcessor processor;
  StructurePreparationData data;
  StructurePreparationFiles testFiles;
  std::string systemName = "plastocyanin";
  Core::Log silentLogger = Core::Log::silent();
  std::shared_ptr<Utils::Settings> settings = std::make_shared<StructurePreparationSettings>();

  void SetUp() override {
    processor.setLog(silentLogger);
    testFiles.workingDirectory = Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName);
    testFiles.initialize();

    testFiles.atomicInfoFile =
        Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "atomic_info.dat");

    testFiles.connectivityFile =
        Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "Connectivity.dat");
  }
};

TEST_F(AStructurePreparationTest, FilesCanBeInitializedWithDefaults) {
  StructurePreparationFiles defaultFiles;
  defaultFiles.initialize();
  ASSERT_THAT(defaultFiles.proteinFile, Eq("rmc.pdb"));
  ASSERT_THAT(defaultFiles.nonRegContainerFile, Eq("nonregular_container.xyz"));
  // now with a prepending working directory
  std::string workDir = Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName);
  defaultFiles.workingDirectory = workDir;
  defaultFiles.initialize();
  ASSERT_THAT(defaultFiles.proteinFile, Eq(Utils::NativeFilenames::combinePathSegments(workDir, "rmc.pdb")));
}

TEST_F(AStructurePreparationTest, ProteinAndNonRegContainerAreSeparatedCorrectly) {
  testFiles.nonRegContainerFile =
      Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "custom_nonRegContainer.xyz");
  testFiles.proteinFile =
      Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "custom_protein.pdb");

  processor.setFiles(testFiles);
  processor.analyzeStructure(plastocyanin_pdb_file);

  const Utils::AtomCollection nonRegContainer = Utils::ChemicalFileHandler::read(testFiles.nonRegContainerFile).first;
  const Utils::PdbStreamHandler handler;
  std::ifstream input;
  input.open(testFiles.proteinFile);
  auto protein = handler.read(input);
  ASSERT_THAT(nonRegContainer.size(), Eq(27));
  ASSERT_THAT(protein[0].size(), Eq(711));

  boost::filesystem::remove_all(testFiles.nonRegContainerFile);
  boost::filesystem::remove_all(testFiles.proteinFile);
}

TEST_F(AStructurePreparationTest, ProteinAndNonRegContainerAreProtonatedCorrectly) {
  // if protein and nonRegContainer file are not present, it should throw an error
  ASSERT_THROW(processor.protonate(), std::runtime_error);
  // set custom files
  testFiles.protonatedProteinFile =
      Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "protonated_protein.xyz");
  testFiles.protonatedNonRegContainerFile =
      Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "protonated_nonRegContainer.xyz");

  processor.setFiles(testFiles);
  if (!Utils::OpenBabelStreamHandler::checkForBinary()) {
    ASSERT_THROW(processor.protonate(), std::runtime_error);
  }
  else {
    processor.protonate();

    ASSERT_TRUE(boost::filesystem::exists(testFiles.protonatedProteinFile));
    ASSERT_TRUE(boost::filesystem::exists(testFiles.protonatedNonRegContainerFile));
    auto protonatedSubstructure = Utils::ChemicalFileHandler::read(testFiles.protonatedProteinFile).first;
    auto protonatedNonRegContainerSubstructure =
        Utils::ChemicalFileHandler::read(testFiles.protonatedNonRegContainerFile).first;
    ASSERT_THAT(protonatedSubstructure.size(), Eq(1414));
    ASSERT_THAT(protonatedNonRegContainerSubstructure.size(), Eq(49));

    boost::filesystem::remove_all(testFiles.protonatedNonRegContainerFile);
    boost::filesystem::remove_all(testFiles.protonatedProteinFile);
  }
}

TEST_F(AStructurePreparationTest, FinalStructureGenerationWorks) {
  ASSERT_THROW(processor.finalize(), std::runtime_error);
  testFiles.systemFile =
      Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "custom_system.xyz");
  testFiles.titrationSitesFile =
      Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "custom_titration_sites.xyz");
  processor.setFiles(testFiles);
  processor.finalize();
  Utils::AtomCollection totalStructure = Utils::ChemicalFileHandler::read(testFiles.systemFile).first;

  // check that files have been written correctly
  ASSERT_TRUE(boost::filesystem::exists(testFiles.systemFile));
  ASSERT_TRUE(boost::filesystem::exists(testFiles.connectivityFile));

  ASSERT_THAT(totalStructure.size(), Eq(1448));

  // check that solvation can be added
  StructureProcessor processor2;
  processor2.settings().modifyBool(SwooseUtilities::SettingsNames::solvateStructure, true);
  processor2.settings().modifyBool(SwooseUtilities::SettingsNames::solvateStructure, 1);
  processor2.setFiles(testFiles);
  processor2.setLog(silentLogger);
  processor2.finalize();

  Utils::AtomCollection totalStructureSolvated = Utils::ChemicalFileHandler::read(testFiles.systemFile).first;
  // For plastocyanin, dont check the exact size of the molecule (because this may change with any change in the
  // solvation protocol in utils, so just check that it is larger than the dry molecule
  bool structureIsLarger = totalStructureSolvated.size() > totalStructure.size();
  ASSERT_TRUE(structureIsLarger);

  StructureProcessor processor3;
  processor3.setLog(silentLogger);
  // Check if changes in the nonRegContainer are accepted correctly
  testFiles.protonatedNonRegContainerFile = Utils::NativeFilenames::combinePathSegments(
      structure_preparation_dir, systemName, "nonregular_container_H_usermodified.xyz");
  processor3.setFiles(testFiles);
  processor3.finalize();
  Utils::AtomCollection totalStructureDifferent = Utils::ChemicalFileHandler::read(testFiles.systemFile).first;
  ASSERT_THAT(totalStructureDifferent.size(), Eq(1450));

  // check that the atomic info file is empty, because all sites are uncharged
  ASSERT_TRUE(boost::filesystem::is_empty(testFiles.atomicInfoFile));

  // check that a file with titrable sites is generated
  ASSERT_TRUE(boost::filesystem::exists(testFiles.titrationSitesFile));

  boost::filesystem::remove_all(testFiles.systemFile);
  boost::filesystem::remove_all(testFiles.connectivityFile);
  boost::filesystem::remove_all(testFiles.atomicInfoFile);
  // boost::filesystem::remove_all(testFiles.titrationSitesFile);
}

TEST_F(AStructurePreparationTest, CustomInformationForNonRegContainerCanBeAdded) {
  // Write custom atomic info file
  std::string customFile = "atomic_info_copper.dat";
  std::ofstream customInfoFile(customFile);
  testFiles.nonRegContainerInfoFile = customFile;
  customInfoFile << "0 2 0\n";
  customInfoFile.close();

  auto fullStructure = Utils::ChemicalFileHandler::read(
      Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "system.xyz"));

  StructurePreparationHelper::mergeProteinAndNonRegContainer(data, testFiles);
  StructurePreparationHelper::handleBoundariesBetweenProteinAndNonRegContainer(data);
  // write the atomic info file for the protein first
  StructurePreparationIO::writeAtomicInfoFileForProtein(data, testFiles.atomicInfoFile);
  StructurePreparationIO::addAtomicInformationForNonRegContainer(testFiles, data.subsystemMapping);

  // StructurePreparationIO::addAtomicInformationForNonRegContainer(testFiles, helper.getSubsystemMapping());

  ASSERT_TRUE(boost::filesystem::exists(testFiles.atomicInfoFile));
  std::ifstream infoFile(testFiles.atomicInfoFile);
  std::string line;
  // the copper atom is atom nr. 1402 in the full system
  std::getline(infoFile, line);
  ASSERT_THAT(line, Eq("1402 2 0"));
  boost::filesystem::remove_all(testFiles.atomicInfoFile);
}

TEST_F(AStructurePreparationTest, TitrableSitesAreCorrectlyDetected) {
  auto structure = Utils::ChemicalFileHandler::read(plastocyanin_H_xyz_file).first;
  StructurePreparationHelper::updatePdbPreparationData(data, structure);
  // Fill the data object
  data.vectorOfNonRegContainerIndices = {1391, 1392, 1393, 1394, 1395, 1396, 1397, 1398, 1399, 1400,
                                         1401, 1402, 1403, 1404, 1405, 1406, 1407, 1408, 1409, 1410,
                                         1411, 1412, 1413, 1414, 1415, 1416, 1417, 1418, 1419, 1420,
                                         1421, 1422, 1423, 1424, 1425, 1426, 1427, 1428, 1429, 1430};
  auto nTitrableSites = StructurePreparationHelper::collectTitrableSites(data);

  ASSERT_THAT(nTitrableSites.size(), Eq(25));
  int numGlu = std::count_if(nTitrableSites.begin(), nTitrableSites.end(),
                             [&](const TitrableSite& s) { return s.residueName == "GLU"; });
  int numAsp = std::count_if(nTitrableSites.begin(), nTitrableSites.end(),
                             [&](const TitrableSite& s) { return s.residueName == "ASP"; });
  int numArg = std::count_if(nTitrableSites.begin(), nTitrableSites.end(),
                             [&](const TitrableSite& s) { return s.residueName == "ARG"; });
  int numLys = std::count_if(nTitrableSites.begin(), nTitrableSites.end(),
                             [&](const TitrableSite& s) { return s.residueName == "LYS"; });
  int numHis = std::count_if(nTitrableSites.begin(), nTitrableSites.end(),
                             [&](const TitrableSite& s) { return s.residueName == "HIS"; });
  int numCys = std::count_if(nTitrableSites.begin(), nTitrableSites.end(),
                             [&](const TitrableSite& s) { return s.residueName == "CYS"; });
  int numTyr = std::count_if(nTitrableSites.begin(), nTitrableSites.end(),
                             [&](const TitrableSite& s) { return s.residueName == "TYR"; });
  ASSERT_THAT(numArg, Eq(0));
  ASSERT_THAT(numCys, Eq(0));
  ASSERT_THAT(numHis, Eq(0));
  ASSERT_THAT(numTyr, Eq(3));
  ASSERT_THAT(numLys, Eq(6));
  ASSERT_THAT(numAsp, Eq(7));
  ASSERT_THAT(numGlu, Eq(9));

  for (auto& site : nTitrableSites) {
    ASSERT_THAT(site.atoms.size(), Eq(site.indicesInFullStructure.size()));
  }

  ASSERT_TRUE(nTitrableSites[0].residueName == "GLU");
  ASSERT_TRUE(nTitrableSites[1].residueName == "ASP");
  ASSERT_TRUE(nTitrableSites[6].residueName == "LYS");
  ASSERT_TRUE(nTitrableSites[17].residueName == "TYR");

  ASSERT_TRUE(nTitrableSites[0].isAcid);
  ASSERT_TRUE(nTitrableSites[1].isAcid);
  ASSERT_FALSE(nTitrableSites[6].isAcid);
  ASSERT_TRUE(nTitrableSites[17].isAcid);

  ASSERT_FALSE(nTitrableSites[0].isBase);
  ASSERT_FALSE(nTitrableSites[1].isBase);
  ASSERT_TRUE(nTitrableSites[6].isBase);
  ASSERT_FALSE(nTitrableSites[17].isBase);

  std::vector<int> correctIndicesInFullStructure = {7, 8, 9, 10, 11, 12, 13, 14, 15, 719, 720, 721, 722, 723};
  ASSERT_THAT(correctIndicesInFullStructure.size(), Eq(nTitrableSites[0].indicesInFullStructure.size()));
  ASSERT_TRUE(std::equal(correctIndicesInFullStructure.begin(), correctIndicesInFullStructure.end(),
                         nTitrableSites[0].indicesInFullStructure.begin()));

  // check on some examples that critical atoms are detected correctly
  // GLU
  ASSERT_THAT(nTitrableSites[0].criticalAtom, Eq(15));
  // ASP
  ASSERT_THAT(nTitrableSites[1].criticalAtom, Eq(54));
  // LYS
  ASSERT_THAT(nTitrableSites[6].criticalAtom, Eq(215));
  // TYR
  ASSERT_THAT(nTitrableSites[17].criticalAtom, Eq(506));
}

TEST_F(AStructurePreparationTest, AminoAcidsAreProtonatedCorrectly) {
  Utils::PdbStreamHandler pdbHandler;
  std::vector<std::string> listOfExclusions = {"PYL", "SEC", "DCY", "MSE"};
  testFiles.protonatedProteinFile =
      Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "protonated_protein.xyz");
  for (const auto& aa : Scine::StructurePreparation::AminoAcids::aminoAcidHierarchy) {
    if (std::find(listOfExclusions.begin(), listOfExclusions.end(), aa) == listOfExclusions.end()) {
      auto aminoAcidInfo = Scine::Swoose::AminoAcidDataForTests::getAminoAcid(aa);
      std::stringstream in(aminoAcidInfo.first);
      auto structure = pdbHandler.read(in);
      ASSERT_THAT(structure.size(), Eq(1));

      ProtonationHandler handler(data, testFiles, settings);
      handler.setProteinStructure(structure[0]);

      // protonate the amino acid
      handler.protonateAllAminoAcids();

      ASSERT_TRUE(boost::filesystem::exists(testFiles.protonatedProteinFile));
      auto protonatedStructure = Utils::ChemicalFileHandler::read(testFiles.protonatedProteinFile).first;
      // check if the correct number of protons is added
      ASSERT_THAT(protonatedStructure.size(), Eq(aminoAcidInfo.second));
      boost::filesystem::remove_all(testFiles.protonatedProteinFile);
    }
  }
}

} // namespace Tests
} // namespace Scine
