/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TestUtilities/ReferenceDataForTests.h"
#include <Core/Log.h>
#include <Database/Collection.h>
#include <Database/Manager.h>
#include <Database/Objects/Calculation.h>
#include <Database/Objects/DenseMatrixProperty.h>
#include <Database/Objects/Model.h>
#include <Database/Objects/NumberProperty.h>
#include <Database/Objects/SparseMatrixProperty.h>
#include <Database/Objects/Structure.h>
#include <Database/Objects/VectorProperty.h>
#include <Swoose/MMParametrization/MMParametrizationSettings.h>
#include <Swoose/MMParametrization/MolecularSystemPartitioner.h>
#include <Swoose/MMParametrization/ParametrizationData.h>
#include <Swoose/MMParametrization/ParametrizationUtils/ConnectivityGenerator.h>
#include <Swoose/MMParametrization/Parametrizer.h>
#include <Swoose/MMParametrization/ReferenceCalculationHelpers/DatabaseHelper.h>
#include <Swoose/QMMM/QmRegionSelection/QmRegionSelector.h>
#include <Swoose/QMMM/QmRegionSelection/QmRegionSelectorSettings.h>
#include <Swoose/QMMM/QmRegionSelection/QmmmReferenceDataManager.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondDetector.h>
#include <Utils/ExternalQC/Orca/OrcaHessianOutputParser.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <gmock/gmock.h>
#include <boost/filesystem.hpp>
#include <thread>

using namespace testing;
namespace Scine {
using namespace MMParametrization;
namespace Tests {

class DatabaseParametrizationTests : public Test {
 public:
  Database::Manager db;
  std::string hostname;
  const std::string testDatabaseName = "test_database_mode";

  Parametrizer parametrizer;
  ParametrizationData data;
  Qmmm::QmmmData qmmmData;
  Core::Log silentLogger = Core::Log::silent();
  const std::string parFile = "Parameters.dat";
  const std::string connFile = "Connectivity.dat";

  void SetUp() override {
    hostname = getIp();
    Database::Credentials credentials(hostname, 27017, testDatabaseName);
    db.setCredentials(credentials);
    srand(42);
  }

  void TearDown() override {
    boost::filesystem::remove(parFile);
    boost::filesystem::remove(connFile);
  }

  std::string getIp() {
    std::string ip;
    if (std::getenv("TEST_MONGO_DB_IP") != nullptr) {
      ip = std::getenv("TEST_MONGO_DB_IP");
    }
    else
      ip = "127.0.0.1";
    return ip;
  }

  bool databaseIsAvailable() {
    Core::Log log;
    std::string ip = getIp();
    Database::Credentials credentials(ip, 27017, testDatabaseName);
    Database::Manager db;
    db.setCredentials(credentials);

    try {
      db.connect();
      return db.isConnected() ? true : false;
    }
    catch (const std::exception& e) {
      log.warning << "Database is not available. Skipping this test... " << Core::Log::endl;
      return false;
    }
  }
};

TEST_F(DatabaseParametrizationTests, ParametrizationFromDatabaseWorks) {
  if (!databaseIsAvailable())
    GTEST_SKIP();

  db.connect();
  db.wipe();
  db.init();

  Utils::AtomCollection smallStructure = Utils::ChemicalFileHandler::read(alanine_xyz_file).first;
  auto settings = std::make_shared<MMParametrizationSettings>();

  data.fullStructure = std::move(smallStructure);
  data.numberOfAtoms = data.fullStructure.size();
  data.vectorOfHessians.resize(1);

  ConnectivityGenerator connectivityGenerator(data, settings, silentLogger);
  connectivityGenerator.generateInitialListsOfNeighbors();

  settings->modifyString(SwooseUtilities::SettingsNames::referenceDataMode, SwooseUtilities::OptionNames::databaseMode);
  settings->modifyString(SwooseUtilities::SettingsNames::databaseHost, hostname);
  settings->modifyInt(SwooseUtilities::SettingsNames::databasePort, 27017);
  settings->modifyString(SwooseUtilities::SettingsNames::databaseName, testDatabaseName);
  settings->modifyInt(SwooseUtilities::SettingsNames::numberAtomsThreshold, 120);
  settings->modifyInt(SwooseUtilities::SettingsNames::databaseSleepTime, 1);

  auto properties = db.getCollection("properties");
  auto calculations = db.getCollection("calculations");
  auto structures = db.getCollection("structures");

  std::string hessianQuery = "{ \"job.order\" : \"scine_hessian\" }";
  std::string structureOptimizationQuery = " { \"job.order\" : \"orca_geometry_optimization\" }";
  std::string bondOrdersQuery = " { \"job.order\" : \"scine_bond_orders\" }";

  MolecularSystemPartitioner partitioner(data, settings, silentLogger);
  partitioner.divideSystem();

  DatabaseHelper helper(data, settings, silentLogger);

  std::vector<Database::Property> storedProperties;
  std::vector<Database::Calculation> vectorOfSubmittedCalculations;

  // Set the optimized structure manually
  Utils::AtomCollection atoms =
      Utils::ChemicalFileHandler::read(Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "0", "opt.xyz"))
          .first;

  // Set the charges manually
  Eigen::VectorXd chargesVector(13);
  chargesVector << 0.099864, -0.227143, 0.092378, 0.032312, 0.030376, -0.077495, 0.030496, 0.024756, 0.029539, 0.214499,
      -0.265597, -0.163429, 0.180475;

  // Set the bond orders manually
  Eigen::SparseMatrix<double> boMatrix = Utils::BondDetector::detectBonds(atoms).getMatrix();

  // Set the hessian matrix manually
  Eigen::MatrixXd hessianMatrix = Utils::ExternalQC::OrcaHessianOutputParser::getHessian(
      Utils::NativeFilenames::combinePathSegments(alanine_ref_calc_dir, "0", "orca.hess"));

  auto executeParametrization = [&]() { helper.runCalculationsAndCollectResults(); };

  auto fillDatabaseWithResults = [&]() {
    while ((calculations->count(structureOptimizationQuery) == 0) && (calculations->count(bondOrdersQuery) == 0))
      std::this_thread::sleep_for(std::chrono::milliseconds(500));

    for (auto iter = calculations->iteratorQuery<Database::Calculation>("{}"); !iter.done(); iter++) {
      auto calc = *iter;
      vectorOfSubmittedCalculations.push_back(calc);
      calc.link(calculations);
      auto structure = Database::Structure(calc.getStructures().at(0));
      structure.link(structures);
      auto model = calc.getModel();
      auto job = calc.getJob();

      auto newLabel = Database::Structure::LABEL::MINIMUM_OPTIMIZED;

      calc.setModel(model);

      auto results = calc.getResults();

      Database::Structure newStructure;
      newStructure.link(structures);

      Database::NumberProperty energy;
      energy.link(properties);
      energy.create(model, "electronic_energy", -323.142716006307);
      storedProperties.push_back(energy);

      if (job.order == "orca_geometry_optimization") {
        Database::VectorProperty hirshfeldCharges;
        hirshfeldCharges.link(properties);
        hirshfeldCharges.create(model, "hirshfeld_charges", chargesVector);
        storedProperties.push_back(hirshfeldCharges);

        newStructure.create(atoms, 0, 1, model, newLabel);
        newStructure.setComment("0");
        newStructure.addProperty("electronic_energy", energy.id());
        newStructure.addProperty("hirshfeld_charges", hirshfeldCharges.id());
        results.properties.push_back(energy.id());
        results.properties.push_back(hirshfeldCharges.id());

        results.structures.push_back(newStructure.id());
        calc.setResults(results);
        calc.setStatus(Database::Calculation::STATUS::COMPLETE);
      }

      else if (job.order == "scine_bond_orders") {
        Database::SparseMatrixProperty bondOrders;
        bondOrders.link(properties);
        bondOrders.create(model, "bond_orders", boMatrix);
        storedProperties.push_back(bondOrders);

        structure.addProperty("electronic_energy", energy.id());
        structure.addProperty("bond_orders", bondOrders.id());
        results.properties.push_back(energy.id());
        results.properties.push_back(bondOrders.id());
        calc.setResults(results);
        calc.setStatus(Database::Calculation::STATUS::COMPLETE);
      }
    }
    while (calculations->count(hessianQuery) == 0)
      std::this_thread::sleep_for(std::chrono::milliseconds(500));

    auto hessianVec = calculations->query<Database::Calculation>(hessianQuery);
    auto calc = hessianVec.at(0);
    vectorOfSubmittedCalculations.push_back(calc);
    calc.link(calculations);
    calc.setModel(calc.getModel());

    auto results = calc.getResults();

    Database::DenseMatrixProperty hessian;
    hessian.link(properties);
    hessian.create(calc.getModel(), "hessian", hessianMatrix);
    storedProperties.push_back(hessian);
    auto structure = Database::Structure(calc.getStructures().at(0));
    structure.link(structures);
    structure.addProperty("hessian", hessian.id());
    results.properties.push_back(hessian.id());
    calc.setResults(results);
    calc.setStatus(Database::Calculation::STATUS::COMPLETE);
  };

  // Start executive thread
  std::thread parametrizationThread(executeParametrization);
  std::thread fillDatabaseThread(fillDatabaseWithResults);

  parametrizationThread.join();
  fillDatabaseThread.join();

  for (auto property : storedProperties)
    ASSERT_TRUE(property.hasId());

  for (auto calculation : vectorOfSubmittedCalculations) {
    calculation.link(calculations);
    ASSERT_EQ(calculation.getStatus(), Database::Calculation::STATUS::ANALYZED);
  }

  auto numberOfStructures = structures->count("{}");
  auto numberOfAnalyzedCalculations = calculations->count("{ \"status\" : \"analyzed\" }");
  auto numberOfStructOptCalculations = calculations->count(structureOptimizationQuery);
  auto numberOfBondOrdersCalculations = calculations->count(bondOrdersQuery);
  auto numberOfHessianCalculations = calculations->count(hessianQuery);
  ASSERT_EQ(numberOfStructures, 2);
  ASSERT_EQ(numberOfStructOptCalculations, 1);
  ASSERT_EQ(numberOfBondOrdersCalculations, 1);
  ASSERT_EQ(numberOfHessianCalculations, 1);
  ASSERT_EQ(numberOfAnalyzedCalculations, 3);

  ASSERT_EQ(data.atomicChargesForEachFragment.at(0).size(), chargesVector.size());

  auto dataStructure = *data.vectorOfOptimizedStructures.at(0);
  ASSERT_EQ(dataStructure.size(), atoms.size());
  Utils::PositionCollection p1 = dataStructure.getPositions();
  Utils::PositionCollection p2 = atoms.getPositions();
  for (int i = 0; i < p1.rows(); ++i) {
    ASSERT_TRUE(Utils::ElementInfo::Z(dataStructure.getElement(i)) == Utils::ElementInfo::Z(atoms.getElement(i)));
    for (int j = 0; j < p1.cols(); ++j) {
      ASSERT_THAT(p1(i, j), DoubleNear(p2(i, j), 1e-4));
    }
  }

  auto dataHessian = *data.vectorOfHessians.at(0);
  ASSERT_THAT(dataHessian.rows(), Eq(hessianMatrix.rows()));
  ASSERT_THAT(dataHessian.cols(), Eq(hessianMatrix.cols()));

  // Start parametrization with reuse_database option
  parametrizer.setLog(silentLogger);
  parametrizer.settings() = *settings;
  parametrizer.settings().modifyBool(SwooseUtilities::SettingsNames::reuseDatabaseKey, true);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, parFile);
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, connFile);

  parametrizer.parametrize(data.fullStructure);

  ASSERT_TRUE(boost::filesystem::exists(parFile));
  ASSERT_TRUE(boost::filesystem::exists(connFile));

  helper.dropDatabase();
  auto numberOfCalculations = calculations->count("{}");
  ASSERT_EQ(numberOfCalculations, 0);
}

TEST_F(DatabaseParametrizationTests, QmRegionSelectionInDatabaseModeWorks) {
  if (!databaseIsAvailable())
    GTEST_SKIP();

  db.connect();
  db.wipe();
  db.init();

  Qmmm::QmRegionSelector qmRegionSelector;
  qmRegionSelector.setLog(silentLogger);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCenterAtom, 302);
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::initialRadiusForQmRegionSelection, 2.0);
  qmRegionSelector.settings().modifyDouble(SwooseUtilities::SettingsNames::cuttingProbability, 0.99);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMinSize, 20);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionCandidateMaxSize, 40);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::qmRegionRefMaxSize, 60);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::numAttemptsPerRadius, 3);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::maxNumRefModels, 1);
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::connectivityFilePath, insulin_connectivity_file);
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::referenceDataMode,
                                           SwooseUtilities::OptionNames::databaseMode);
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::databaseHost, hostname);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::databasePort, 27017);
  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::databaseName, testDatabaseName);
  qmRegionSelector.settings().modifyInt(SwooseUtilities::SettingsNames::databaseSleepTime, 1);

  // Read in structure
  Utils::AtomCollection structure = Utils::ChemicalFileHandler::read(insulin_xyz_file).first;

  // Should raise an error when SFAM parameter file is not available
  ASSERT_THROW(qmRegionSelector.generateQmRegion(structure), std::runtime_error);

  qmRegionSelector.settings().modifyString(SwooseUtilities::SettingsNames::parameterFilePath, insulin_parameter_file);

  auto properties = db.getCollection("properties");
  auto calculations = db.getCollection("calculations");
  auto structures = db.getCollection("structures");

  std::string qmmmQuery = "{ \"job.order\" : \"swoose_qmmm_forces\" }";

  auto executeQmRegionSelection = [&]() { qmRegionSelector.generateQmRegion(structure); };

  int numberOfNewCalculations = 0;
  int numberOfHoldCalculations = 0;
  int sizeOfSelectedQmRegionIndices = 0; // will later be set equal to the correct value

  auto fillDatabaseWithResults = [&]() {
    while ((calculations->count(qmmmQuery) == 0))
      std::this_thread::sleep_for(std::chrono::milliseconds(500));

    std::this_thread::sleep_for(std::chrono::milliseconds(500)); // make sure all calculations are in the db already

    std::vector<Swoose::ForcesCollection> refForces =
        Swoose::ReferenceDataForTests::parseReferenceForces(reference_forces_insulin_file, 2, 328);

    int lastIndex = calculations->count(qmmmQuery) - 1;
    std::string lastIndexAsString = std::to_string(lastIndex);

    bool selectionWasAnalyzed = false; // tracks whether we already processed the later selected QM region
    for (auto iter = calculations->iteratorQuery<Database::Calculation>(qmmmQuery); !iter.done(); iter++) {
      auto calc = *iter;
      calc.link(calculations);
      auto structure = Database::Structure(calc.getStructures()[0]);
      structure.link(structures);
      auto model = calc.getModel();
      auto job = calc.getJob();

      calc.setModel(model);

      // Get all the qm atoms for this structure
      std::vector<int> listOfQmAtoms = calc.getSetting(SwooseUtilities::SettingsNames::qmAtomsList);

      auto results = calc.getResults();
      if (calc.getStatus() == Database::Calculation::STATUS::HOLD) {
        numberOfHoldCalculations++;
      }
      else if (calc.getStatus() == Database::Calculation::STATUS::NEW) {
        numberOfNewCalculations++;
        if (!selectionWasAnalyzed || calc.getComment() == lastIndexAsString) {
          if (!selectionWasAnalyzed) {
            if (calc.getComment() != lastIndexAsString) { // the selected candidate should not be the reference model
              sizeOfSelectedQmRegionIndices = listOfQmAtoms.size();
              selectionWasAnalyzed = true;
            }
          }

          Database::NumberProperty energy;
          energy.link(properties);
          energy.create(model, "electronic_energy", -217.5043024);
          energy.setCalculation(calc.id());
          energy.setStructure(structure.id());

          Database::DenseMatrixProperty atomicForces;
          atomicForces.link(properties);
          atomicForces.create(model, "atomic_forces", refForces.at(0));
          atomicForces.setCalculation(calc.id());
          atomicForces.setStructure(structure.id());

          structure.addProperty("atomic_forces", atomicForces.id());
          structure.addProperty("electronic_energy", energy.id());
          results.properties.push_back(energy.id());
          results.properties.push_back(atomicForces.id());
          calc.setResults(results);

          calc.setStatus(Database::Calculation::STATUS::COMPLETE);
        }
        else {
          Database::NumberProperty energy2;
          energy2.link(properties);
          energy2.create(model, "electronic_energy", -261.858301);
          energy2.setCalculation(calc.id());
          energy2.setStructure(structure.id());

          Database::DenseMatrixProperty atomicForces2;
          atomicForces2.link(properties);
          atomicForces2.create(model, "atomic_forces", refForces.at(1));
          atomicForces2.setCalculation(calc.id());
          atomicForces2.setStructure(structure.id());

          structure.addProperty("atomic_forces", atomicForces2.id());
          structure.addProperty("electronic_energy", energy2.id());
          results.properties.push_back(energy2.id());
          results.properties.push_back(atomicForces2.id());
          calc.setResults(results);

          calc.setStatus(Database::Calculation::STATUS::COMPLETE);
        }
      }
    }
  };

  // Start executive thread
  std::thread execute(executeQmRegionSelection);
  // Start thread that fills database with results
  std::thread fillDB(fillDatabaseWithResults);

  execute.join();
  fillDB.join();

  for (auto iter = calculations->iteratorQuery<Database::Calculation>(qmmmQuery); !iter.done(); iter++) {
    auto calc = *iter;
    auto settings = calc.getSettings();
    // check that qm region selection settings are not in the DB
    ASSERT_FALSE(settings.valueExists("qm_region_center_atoms"));
    ASSERT_FALSE(settings.valueExists("initial_radius"));
    ASSERT_FALSE(settings.valueExists("cutting_probability"));
    ASSERT_FALSE(settings.valueExists("qm_region_min_size"));
    ASSERT_FALSE(settings.valueExists("qm_region_max_size"));
    ASSERT_FALSE(settings.valueExists("ref_max_size"));
    ASSERT_FALSE(settings.valueExists("num_attempts_per_radius"));
    ASSERT_FALSE(settings.valueExists("max_num_ref_models"));
    // check that QM/MM calculation settings are in the DB
    ASSERT_TRUE(settings.valueExists("qm_atoms"));
    ASSERT_TRUE(settings.valueExists("molecular_charge"));
    ASSERT_TRUE(settings.valueExists("spin_multiplicity"));
  }

  auto numberOfStructures = structures->count("{}");
  auto numberOfAnalyzedCalculations = calculations->count("{ \"status\" : \"analyzed\" }");
  auto numberOfQmmmCalculations = calculations->count(qmmmQuery);
  auto numberOfNotPerformedCalculations = calculations->count("{ \"status\" : \"hold\" }");
  ASSERT_EQ(numberOfStructures, numberOfNewCalculations + numberOfHoldCalculations);
  ASSERT_EQ(numberOfAnalyzedCalculations, numberOfNewCalculations);
  ASSERT_EQ(numberOfQmmmCalculations, numberOfNewCalculations + numberOfHoldCalculations);
  ASSERT_EQ(numberOfNotPerformedCalculations, numberOfHoldCalculations);

  auto qmRegion = qmRegionSelector.getQmRegionStructure();
  auto qmRegionIndices = qmRegionSelector.getQmRegionIndices();
  auto qmRegionInfo = qmRegionSelector.getQmRegionChargeAndMultiplicity();

  ASSERT_TRUE(std::abs(qmRegion.size() - sizeOfSelectedQmRegionIndices) < 6); // less than 6 saturating atoms
  ASSERT_TRUE(qmRegion.size() >= sizeOfSelectedQmRegionIndices);
  ASSERT_THAT(qmRegionIndices.size(), Eq(qmRegion.size()));
  ASSERT_THAT(qmRegionInfo.first, Eq(0));
  ASSERT_THAT(qmRegionInfo.second, Eq(1));
}

} // namespace Tests
} // namespace Scine
