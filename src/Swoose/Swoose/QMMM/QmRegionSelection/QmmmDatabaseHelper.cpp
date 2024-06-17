/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
                                       const std::vector<QmmmModel>& qmmmReferenceModels, const QmmmData& qmmmData,
                                       const Utils::Settings& calculatorSettings)
  : log_(log),
    settings_(settings),
    database_(std::make_unique<Database::Manager>()),
    qmmmModelCandidates_(qmmmModelCandidates),
    qmmmReferenceModels_(qmmmReferenceModels),
    structure_(structure),
    bondOrders_(bondOrders),
    qmmmData_(qmmmData),
    calculatorSettings_(calculatorSettings) {
  sleepTime_ = settings.getInt(SwooseUtilities::SettingsNames::databaseSleepTime);
  auto databaseName = settings_.getString(SwooseUtilities::SettingsNames::databaseName);
  auto databasePort = settings_.getInt(SwooseUtilities::SettingsNames::databasePort);
  auto databaseHost = settings_.getString(SwooseUtilities::SettingsNames::databaseHost);
  auto reuseDatabase = settings_.getBool(SwooseUtilities::SettingsNames::reuseDatabaseKey);
  Database::Credentials credentials(databaseHost, databasePort, databaseName);
  database_->setCredentials(credentials);
  try {
    log_.output << "Connecting to database..." << Core::Log::nl << Core::Log::endl;
    database_->connect();
    if (!reuseDatabase) {
      database_->wipe();
    }
    database_->init();
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
  int numCandidateModels = qmmmModelCandidates_.size();
  long unsigned int numModels = numCandidateModels + qmmmReferenceModels_.size();
  std::vector<Database::ID> structureIDs(numModels);
  for (long unsigned int i = 0; i < numModels; ++i) {
    int molecularCharge, spinMultiplicity;
    if (int(i) < numCandidateModels) {
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
    if (settings_.getIntList(SwooseUtilities::SettingsNames::qmRegionCenterAtoms).size() == 1) {
      bool submitted = submitCalculation(candidateModel, structureIDs[modelCounter], modelCounter, collCalc,
                                         qmmmData_.symmetryScores.at(modelCounter), maxAllowedSymmetryScore);
      if (!submitted)
        holdCounter++;
    }
    else
      submitCalculation(candidateModel, structureIDs[modelCounter], modelCounter, collCalc);
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
  auto parameterFilePath = settings_.getString(Utils::SettingsNames::parameterFilePath);
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
  if (yamlNodeQmmmSettings[Utils::SettingsNames::externalProgramNProcs])
    job.cores = yamlNodeQmmmSettings[Utils::SettingsNames::externalProgramNProcs].as<int>();
  else if (calculatorSettings_.valueExists(Utils::SettingsNames::externalProgramNProcs))
    job.cores = calculatorSettings_.getInt(Utils::SettingsNames::externalProgramNProcs);
  else
    job.cores = 1;
  if (yamlNodeQmmmSettings[Utils::SettingsNames::externalProgramMemory])
    job.memory = yamlNodeQmmmSettings[Utils::SettingsNames::externalProgramMemory].as<int>();
  else if (calculatorSettings_.valueExists(Utils::SettingsNames::externalProgramMemory))
    job.memory = calculatorSettings_.getInt(Utils::SettingsNames::externalProgramMemory);
  else
    job.memory = 2.0 + (job.cores / 2) + (qmmmModel.structure.size() / 50.0);
  job.disk = 1.0 + (qmmmModel.structure.size() / 1024.0);

  // handle model fields, not complete yet
  // method family and program are still a workaround around the proper QM/MM search of Puffin
  std::string method = calculatorSettings_.valueExists(Utils::SettingsNames::method)
                           ? calculatorSettings_.getString(Utils::SettingsNames::method)
                           : "any";
  std::string basisSet;
  // if we have not received calculator settings, we guess 'any'
  // if we have received calculator settings, and they do not contain a basis set, we set it to an empty string
  if (calculatorSettings_.empty()) {
    basisSet = "any";
  }
  else {
    basisSet = calculatorSettings_.valueExists(Utils::SettingsNames::basisSet)
                   ? calculatorSettings_.getString(Utils::SettingsNames::basisSet)
                   : "";
  }
  std::string actualProgram = settings_.valueExists(Utils::SettingsNames::program)
                                  ? settings_.getString(Utils::SettingsNames::program)
                                  : "Any/Swoose";
  std::string actualMethodFamily = settings_.valueExists(Utils::SettingsNames::methodFamily)
                                       ? settings_.getString(Utils::SettingsNames::methodFamily)
                                       : "";

  // Remove QM region selection settings from yaml node and give remaining node to the selector
  std::vector<std::string> keysToRemove;
  for (YAML::const_iterator it = yamlNodeQmmmSettings.begin(); it != yamlNodeQmmmSettings.end(); ++it) {
    auto key = it->first.as<std::string>();
    QmRegionSelector qmRegionSelector;
    if (qmRegionSelector.settings().valueExists(key)) {
      keysToRemove.push_back(key);
    }
    else if (key == Utils::SettingsNames::method) {
      method = it->second.as<std::string>();
      keysToRemove.push_back(key);
    }
    else if (key == Utils::SettingsNames::basisSet) {
      basisSet = it->second.as<std::string>();
      keysToRemove.push_back(key);
    }
    else if (key == Utils::SettingsNames::methodFamily) {
      actualMethodFamily = it->second.as<std::string>();
    }
    else if (key == Utils::SettingsNames::program) {
      actualProgram = it->second.as<std::string>();
    }
  }
  if (actualMethodFamily.empty() || actualMethodFamily.find('/') == std::string::npos) {
    throw std::runtime_error("Require a method_family such as 'PM6/SFAM'");
  }
  if (actualProgram.find('/') == std::string::npos) {
    throw std::runtime_error(
        "If you specify a program, it must include the QM and MM program such as 'Sparrow/Swoose'");
  }
  for (const auto& key : keysToRemove)
    yamlNodeQmmmSettings.remove(key);

  Database::Model model("QM/MM", method, basisSet);
  model.program = "swoose"; // do not write actualProgram here due to lacking Puffin features
  calc.create(model, job, {structureID});
  calc.setComment(std::to_string(qmmmModelIndex));
  for (const auto& [k, v] : calculatorSettings_.items()) {
    calc.setSetting(k, v);
  }
  // relies on the fact that the Puffin Swoose job picks these up
  // todo remove once Puffin is properly updated
  calc.setSetting("qm_model", actualMethodFamily.substr(0, actualMethodFamily.find('/')));
  calc.setSetting("qm_module", actualProgram.substr(0, actualProgram.find('/')));
  calc.setSetting(Utils::SettingsNames::methodFamily, actualMethodFamily);
  calc.setSetting(Utils::SettingsNames::program, actualProgram);

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
  calc.setSetting(Utils::SettingsNames::qmAtomsList, std::move(validQmAtoms));

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
