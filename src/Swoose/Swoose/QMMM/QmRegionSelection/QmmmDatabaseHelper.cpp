/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmmmDatabaseHelper.h"
#include "QmRegionSelector.h"
#include "QmRegionSelectorSettings.h"
#include "QmmmReferenceDataManager.h"
#include "SymmetryScores.h"
#include <Core/Log.h>
#include <Database/Collection.h>
#include <Database/Manager.h>
#include <Database/Objects/Calculation.h>
#include <Database/Objects/DenseMatrixProperty.h>
#include <Database/Objects/Model.h>
#include <Database/Objects/SparseMatrixProperty.h>
#include <Database/Objects/StringProperty.h>
#include <Database/Objects/Structure.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <thread>

namespace Scine {
namespace Qmmm {

QmmmDatabaseHelper::QmmmDatabaseHelper(const Utils::Settings& settings, Core::Log& log,
                                       const Utils::AtomCollection& structure, const Utils::BondOrderCollection& bondOrders,
                                       const std::vector<QmmmModel>& qmmmModelCandidates,
                                       const std::vector<QmmmModel>& qmmmReferenceModels, const QmmmData& qmmmData)
  : database_(std::make_unique<Database::Manager>()),
    settings_(settings),
    log_(log),
    structure_(structure),
    qmmmModelCandidates_(qmmmModelCandidates),
    qmmmReferenceModels_(qmmmReferenceModels),
    bondOrders_(bondOrders),
    qmmmData_(qmmmData) {
  sleepTime_ = settings.getInt(SwooseUtilities::SettingsNames::databaseSleepTime);
  auto databaseName = settings_.getString(SwooseUtilities::SettingsNames::databaseName);
  auto databasePort = settings_.getInt(SwooseUtilities::SettingsNames::databasePort);
  auto databaseHost = settings_.getString(SwooseUtilities::SettingsNames::databaseHost);
  Database::Credentials credentials(databaseHost, databasePort, databaseName);
  database_->setCredentials(credentials);
  try {
    log_.output << "Connecting to database..." << Core::Log::nl << Core::Log::endl;
    database_->connect();
  }
  catch (const std::exception& e) {
    throw std::runtime_error("Connecting to the specified database failed!");
  }
}

std::vector<ForcesCollection> QmmmDatabaseHelper::calculateForces() {
  // Get the collections from the database
  auto collCalc = database_->getCollection("calculations");
  auto collStruct = database_->getCollection("structures");
  auto collProp = database_->getCollection("properties");

  // Add parameters properties to database
  Database::Model model("", "", "");
  auto parametersAsString = getParametersAsString();
  Database::StringProperty parametersProperty;
  parametersProperty.link(collProp);
  auto parametersID = parametersProperty.create(model, "sfam_parameters", parametersAsString);

  // Add bond orders properties to database
  Database::SparseMatrixProperty bondOrdersProperty;
  bondOrdersProperty.link(collProp);
  auto bondOrdersID = bondOrdersProperty.create(model, "bond_orders", bondOrders_.getMatrix());

  // Add structure to database
  auto numCandidateModels = qmmmModelCandidates_.size();
  auto numModels = numCandidateModels + qmmmReferenceModels_.size();
  std::vector<Database::ID> structureIDs(numModels);
  for (int i = 0; i < numModels; ++i) {
    int molecularCharge, spinMultiplicity;
    if (i < numCandidateModels) {
      molecularCharge = qmmmModelCandidates_.at(i).molecularCharge;
      spinMultiplicity = qmmmModelCandidates_.at(i).spinMultiplicity;
    }
    else {
      molecularCharge = qmmmReferenceModels_.at(i - numCandidateModels).molecularCharge;
      spinMultiplicity = qmmmReferenceModels_.at(i - numCandidateModels).spinMultiplicity;
    }
    Database::Structure structure;
    structure.link(collStruct);
    auto structureID = structure.create(structure_, molecularCharge, spinMultiplicity);
    auto label = Database::Structure::LABEL::MINIMUM_GUESS;
    structure.setLabel(label);
    structure.addProperty("sfam_parameters", parametersID);
    structure.addProperty("bond_orders", bondOrdersID);
    structureIDs[i] = structureID;
  }

  const double maxAllowedSymmetryScore = calculateMaxAllowedSymmetryScore(qmmmData_, settings_);
  int modelCounter = 0;
  int holdCounter = 0;
  for (const auto& candidateModel : qmmmModelCandidates_) {
    bool submitted = submitCalculation(candidateModel, structureIDs[modelCounter], modelCounter, collCalc,
                                       qmmmData_.symmetryScores.at(modelCounter), maxAllowedSymmetryScore);
    if (!submitted)
      holdCounter++;
    modelCounter++;
  }
  for (const auto& refModel : qmmmReferenceModels_) {
    submitCalculation(refModel, structureIDs[modelCounter], modelCounter, collCalc);
    modelCounter++;
  }

  log_.debug << "Number of model candidates that are ignored due to a large symmetry score: " << holdCounter
             << Core::Log::nl << Core::Log::endl;

  bool unfinishedCalculationsExist = true;
  while (unfinishedCalculationsExist) {
    // Verify connection to database
    while (!database_->isConnected()) {
      log_.output << "Cannot connect to database. Trying again in " << sleepTime_ << " seconds." << Core::Log::nl
                  << Core::Log::endl;
      std::this_thread::sleep_for(std::chrono::seconds(sleepTime_));
    }

    const std::string incompleteCalculationsQuery =
        "{ \"$or\" : [{ \"status\" : \"new\" }, { \"status\" : \"pending\" }] }";
    auto incompleteVec = collCalc->query<Database::Calculation>(incompleteCalculationsQuery);
    unfinishedCalculationsExist = !incompleteVec.empty();
    if (unfinishedCalculationsExist) {
      log_.debug << "Sleeping for " << sleepTime_ << " seconds..." << Core::Log::nl << Core::Log::endl;
      std::this_thread::sleep_for(std::chrono::seconds(sleepTime_));
    }
  }

  return collectCalculationResults(collCalc, collStruct, collProp);
}

std::string QmmmDatabaseHelper::getParametersAsString() const {
  auto parameterFilePath = settings_.getString(SwooseUtilities::SettingsNames::parameterFilePath);
  std::ifstream file(parameterFilePath);
  if (!file.is_open())
    throw std::runtime_error("The parameter file " + parameterFilePath + " cannot be opened.");
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

bool QmmmDatabaseHelper::submitCalculation(const QmmmModel& qmmmModel, const Database::ID& structureID,
                                           int qmmmModelIndex, std::shared_ptr<Database::Collection> calculations,
                                           double symmetryScore, double maxAllowedSymmetryScore) {
  // Get yaml settings for QM/MM
  auto yamlSettingsPath = settings_.getString(SwooseUtilities::SettingsNames::yamlSettingsFilePath);
  YAML::Node yamlNodeQmmmSettings;
  if (!yamlSettingsPath.empty())
    yamlNodeQmmmSettings = YAML::LoadFile(yamlSettingsPath);

  Database::Calculation calc;
  calc.link(calculations);

  Database::Calculation::Job job("swoose_qmmm_forces");
  if (yamlNodeQmmmSettings["external_program_nprocs"])
    job.cores = yamlNodeQmmmSettings["external_program_nprocs"].as<int>();
  else
    job.cores = 1;
  if (yamlNodeQmmmSettings["external_program_memory"])
    job.memory = yamlNodeQmmmSettings["external_program_memory"].as<int>();
  else
    job.memory = 2.0 + (job.cores / 2) + (qmmmModel.structure.size() / 50.0);
  job.disk = 2.0 + (qmmmModel.structure.size() / 50.0);

  Database::Model model("QM/MM", "any", "any");
  model.program = "swoose";
  calc.create(model, job, {structureID});
  calc.setComment(std::to_string(qmmmModelIndex));

  // Remove QM region selection settings from yaml node and give remaining node to the selector
  std::vector<std::string> keysToRemove;
  for (YAML::const_iterator it = yamlNodeQmmmSettings.begin(); it != yamlNodeQmmmSettings.end(); ++it) {
    auto key = it->first.as<std::string>();
    QmRegionSelector qmRegionSelector;
    if (qmRegionSelector.settings().valueExists(key)) {
      keysToRemove.push_back(key);
    }
  }
  for (const auto& key : keysToRemove)
    yamlNodeQmmmSettings.remove(key);

  if (yamlNodeQmmmSettings.IsMap()) {
    auto settingsAsValueCollection = Utils::deserializeValueCollection(yamlNodeQmmmSettings);
    for (auto it = settingsAsValueCollection.begin(); it != settingsAsValueCollection.end(); ++it)
      calc.setSetting(it->first, it->second);
  }

  calc.setSetting(Utils::SettingsNames::spinMultiplicity, qmmmModel.spinMultiplicity);
  calc.setSetting(Utils::SettingsNames::molecularCharge, qmmmModel.molecularCharge);

  std::vector<int> validQmAtoms;
  validQmAtoms.reserve(qmmmModel.qmAtomIndices.size());
  std::copy_if(std::begin(qmmmModel.qmAtomIndices), std::end(qmmmModel.qmAtomIndices), std::back_inserter(validQmAtoms),
               [](int a) -> bool { return a >= 0; });
  calc.setSetting(SwooseUtilities::SettingsNames::qmAtomsList, std::move(validQmAtoms));

  if (symmetryScore <= maxAllowedSymmetryScore)
    calc.setStatus(Database::Calculation::STATUS::NEW);
  else {
    calc.setStatus(Database::Calculation::STATUS::HOLD);
    return false;
  }
  return true;
}

std::vector<ForcesCollection>
QmmmDatabaseHelper::collectCalculationResults(std::shared_ptr<Database::Collection> calculations,
                                              std::shared_ptr<Database::Collection> structures,
                                              std::shared_ptr<Database::Collection> properties) {
  std::vector<ForcesCollection> forces(qmmmModelCandidates_.size() + qmmmReferenceModels_.size());
  const std::string queryForCompletedCalculations = "{ \"status\" : \"complete\" }";
  for (auto iter = calculations->iteratorQuery<Database::Calculation>(queryForCompletedCalculations); !iter.done(); iter++) {
    auto calc = *iter;
    calc.link(calculations);
    int modelIndex = std::stoi(calc.getComment());

    // Get structure
    auto structureID = calc.getStructures().at(0);
    auto structure = structures->get<Database::Structure>(structureID);
    structure.link(structures);

    auto forcesID = structure.getProperties("atomic_forces").at(0);
    auto forcesProperty = properties->get<Database::DenseMatrixProperty>(forcesID);
    forcesProperty.link(properties);
    forces.at(modelIndex) = forcesProperty.getData();

    calc.setStatus(Database::Calculation::STATUS::ANALYZED);
  }
  return forces;
}

} // namespace Qmmm
} // namespace Scine
