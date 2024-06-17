/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include <Swoose/MMParametrization/MMParametrizationSettings.h>
#include <Swoose/MMParametrization/ParametrizationData.h>
#include <Swoose/MMParametrization/ParametrizationUtils/TitrationHelper.h>
#include <Swoose/MMParametrization/ReferenceCalculationHelpers/ReferenceCalculationsIO.h>
#include <Swoose/StructurePreparation/StructurePreparationHelper.h>
#include <Swoose/StructurePreparation/StructureProcessor.h>
#include <Swoose/Utilities/TitrationFileHandler.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
using namespace StructurePreparation;
using namespace MMParametrization;
namespace Tests {

/**
 * @class APdbPreparationTest PdbPreparationTest.cpp
 * @class Tests the automated PDB preparation for plastocyanin.
 */

class ATitrationTest : public Test {
 public:
  StructureProcessor preparator;
  ParametrizationData data;
  StructurePreparationFiles files;
  TitrationResults results;
  std::string systemName = "plastocyanin";
  Core::Log silentLogger = Core::Log::silent();
  std::shared_ptr<Utils::Settings> settings = std::make_shared<MMParametrizationSettings>();

  void SetUp() override {
    preparator.setLog(silentLogger);
    files.protonatedProteinFile =
        Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, systemName, "protein_H.pdb");
  }
};

TEST_F(ATitrationTest, TitrationSiteFileIsParsedCorrectly) {
  std::vector<int> listOfCorrectSites = {15,  54,  62,  123, 171, 180, 215, 299, 308, 316, 325, 365, 384,
                                         421, 430, 438, 487, 506, 515, 553, 562, 585, 594, 617, 680};
  std::string filename =
      Utils::NativeFilenames::combinePathSegments(structure_preparation_dir, "plastocyanin", "titrable_sites.dat");
  data.siteIspHSensitive.resize(1447);
  SwooseUtilities::TitrationFileHandler::readTitrationSitesFromFile(results, filename, 1447, 1447,
                                                                    data.pHSensitiveSites, data.siteIspHSensitive);
  ASSERT_THAT(results.sites.size(), Eq(25));
  for (auto& index : listOfCorrectSites) {
    ASSERT_TRUE(data.siteIspHSensitive.at(index));
    ASSERT_THAT(data.pHSensitiveSites.count(index), Eq(1));
  }

  ASSERT_THAT(data.pHSensitiveSites.find(15)->second, Eq("GLU"));
  ASSERT_THAT(data.pHSensitiveSites.find(680)->second, Eq("LYS"));

  ASSERT_THAT(results.sites.at(0).residueName, Eq("GLU"));
  ASSERT_THAT(results.sites.at(24).residueName, Eq("LYS"));

  ASSERT_TRUE(results.sites.at(0).isAcid);
  ASSERT_FALSE(results.sites.at(0).isBase);
}

TEST_F(ATitrationTest, ProtonationStateIsChangesCorrectly) {
  TitrationHelper helper(settings);
  std::string residueName = "GLU";
  bool isBase = false;
  int criticalAtomIndex = 8;
  std::vector<int> superfluousHydrogens = {17};
  auto structure = Utils::ChemicalFileHandler::read(glutamine_xyz_file).first;
  auto deprotStructure =
      helper.changeProtonationState(structure, residueName, isBase, criticalAtomIndex, superfluousHydrogens);
  ASSERT_THAT(structure.size(), Eq(deprotStructure.size() + 1));

  std::stringstream stream("24\n\n"
                           "N    1.620316   -0.364789   3.880711\n"
                           "C    0.698556   -1.010639   2.926220\n"
                           "C    1.513510   -1.592707   1.747226\n"
                           "O    2.797029   -1.192544   1.757088\n"
                           "H    1.258121    0.535294   4.229722\n"
                           "H    2.839127   -0.644849   2.602451\n"
                           "H    0.159512   -1.859267   3.379151\n"
                           "C   -0.327123   -0.008511   2.377357\n"
                           "C   -1.333017    0.476182   3.440012\n"
                           "C   -1.625777    1.983021   3.407511\n"
                           "C   -0.398171    2.885483   3.581186\n"
                           "N    0.427317    2.472301   4.731688\n"
                           "H   -0.845857   -0.477709   1.532231\n"
                           "H    0.235318    0.839942   1.955838\n"
                           "H   -0.958257    0.217591   4.441332\n"
                           "H   -2.279537   -0.071321   3.326108\n"
                           "H   -2.357259    2.216193   4.199299\n"
                           "H   -2.113868    2.257529   2.458910\n"
                           "H   -0.741435    3.934200   3.644647\n"
                           "H    0.243239    2.826116   2.689243\n"
                           "H   -0.115129    2.534089   5.595757\n"
                           "H    1.215641    3.110490   4.844479\n"
                           "O    1.036290   -2.312054   0.895898\n"
                           "H    1.805454   -0.960040   4.685935\n");

  auto baseStructure = Utils::XyzStreamHandler::read(stream);
  Utils::ChemicalFileHandler::write("base.xyz", baseStructure);
  residueName = "LYS";
  isBase = true;
  criticalAtomIndex = 11;
  superfluousHydrogens.clear();
  superfluousHydrogens = {20, 21};
  auto protStructure =
      helper.changeProtonationState(baseStructure, residueName, isBase, criticalAtomIndex, superfluousHydrogens);
  ASSERT_THAT(superfluousHydrogens.size(), Eq(2));
  ASSERT_THAT(baseStructure.size(), Eq(protStructure.size() - 1));
}

TEST_F(ATitrationTest, EnergiesAndStructuresAreParsedCorrectlyInReadMode) {
  int fragmentIndex = 0;
  data.siteIspHSensitive.resize(1);
  data.siteIspHSensitive.at(0) = true;

  ReferenceCalculationsIO::saveAdditionalStructuresForTitration(data, results, fragmentIndex, glutamine_ref_calc_dir);
  auto refStructure = *results.vectorOfOptimizedNonRefStructures.find(0)->second;
  ASSERT_THAT(refStructure.size(), Eq(17));

  bool parseTurbomole = false;
  ReferenceCalculationsIO::parseElectronicEnergiesForTitration(data, results, fragmentIndex, glutamine_ref_calc_dir,
                                                               parseTurbomole, settings);
  auto refEnergy = results.electronicEnergies.find(fragmentIndex)->second.first;
  auto nonRefEnergy = results.electronicEnergies.find(fragmentIndex)->second.second;
  ASSERT_THAT(refEnergy, DoubleNear(-550.078073, 1e-6));
  ASSERT_THAT(nonRefEnergy, DoubleNear(-549.342740, 1e-6));

  settings->modifyBool(SwooseUtilities::SettingsNames::useThermoChemistryForTitration, true);
  ReferenceCalculationsIO::parseElectronicEnergiesForTitration(data, results, fragmentIndex, glutamine_ref_calc_dir,
                                                               parseTurbomole, settings);
  refEnergy = results.electronicEnergies.find(fragmentIndex)->second.first;
  nonRefEnergy = results.electronicEnergies.find(fragmentIndex)->second.second;

  ASSERT_THAT(refEnergy, DoubleNear(-549.97977774, 1e-6));
  ASSERT_THAT(nonRefEnergy, DoubleNear(-549.25720113, 1e-6));

  // works for Turbomole too
  parseTurbomole = true;
  ASSERT_THROW(ReferenceCalculationsIO::parseElectronicEnergiesForTitration(
                   data, results, fragmentIndex, glutamine_ref_calc_dir, parseTurbomole, settings),
               std::runtime_error);
  settings->modifyBool(SwooseUtilities::SettingsNames::useThermoChemistryForTitration, false);
  ReferenceCalculationsIO::parseElectronicEnergiesForTitration(data, results, fragmentIndex, glutamine_ref_calc_dir,
                                                               parseTurbomole, settings);
  refEnergy = results.electronicEnergies.find(fragmentIndex)->second.first;
  nonRefEnergy = results.electronicEnergies.find(fragmentIndex)->second.second;
  ASSERT_THAT(refEnergy, DoubleNear(-550.0778959564, 1e-6));
  ASSERT_THAT(nonRefEnergy, DoubleNear(-549.3425258135, 1e-6));
}

} // namespace Tests
} // namespace Scine
