/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DatabaseHelper.h"
#include "../CalculationManager.h"
#include "../MMParametrizationSettings.h"
#include "../ParametrizationData.h"
#include "../ParametrizationUtils/FragmentDataDistributor.h"
#include "DatabaseJobSubmissionHelper.h"
#include "DatabaseOrderNames.h"
#include <Core/Log.h>
#include <Database/Collection.h>
#include <Database/Manager.h>
#include <Database/Objects/Calculation.h>
#include <Database/Objects/DenseMatrixProperty.h>
#include <Database/Objects/SparseMatrixProperty.h>
#include <Database/Objects/Structure.h>
#include <Database/Objects/VectorProperty.h>
#include <Utils/ExternalQC/Orca/OrcaCalculatorSettings.h>
#include <Utils/Properties/AtomicCharges/ChargeModel5.h>
#include <thread>

namespace Scine {
namespace MMParametrization {

namespace {
static constexpr const char* nameOfBondOrdersProperty = "bond_orders";
static constexpr const char* nameOfHessianProperty = "hessian";
static constexpr const char* nameOfCm5ChargesProperty = "cm5_charges";
static constexpr const char* nameOfHirshfeldChargesProperty = "hirshfeld_charges";
static constexpr const char* nameOfLoewdinChargesProperty = "loewdin_charges";
// The following property name is used if atomic charges are calculated via SCINE single point order:
static constexpr const char* nameOfAtomicChargesProperty = "atomic_charges";
} // namespace

DatabaseHelper::DatabaseHelper(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings, Core::Log& log)
  : data_(data), settings_(settings), database_(std::make_unique<Database::Manager>()), log_(log) {
  numberOfStructures_ = static_cast<int>(data_.vectorOfStructures.size());
  reuseDatabase_ = settings->getBool(SwooseUtilities::SettingsNames::reuseDatabaseKey);
  sleepTime_ = settings->getInt(SwooseUtilities::SettingsNames::databaseSleepTime);
  refineConnectivity_ = settings->getBool(SwooseUtilities::SettingsNames::refineConnectivity);
  useGaussian_ = settings_->getBool(SwooseUtilities::SettingsNames::useGaussianOptionKey);
  earlyTerminationEnabled_ = settings_->getBool(SwooseUtilities::SettingsNames::enableEarlyTerminationKey);
  auto databaseName = settings_->getString(SwooseUtilities::SettingsNames::databaseName);
  auto databasePort = settings_->getInt(SwooseUtilities::SettingsNames::databasePort);
  auto databaseHost = settings_->getString(SwooseUtilities::SettingsNames::databaseHost);
  Database::Credentials credentials(databaseHost, databasePort, databaseName);
  database_->setCredentials(credentials);
  referenceProgram_ = settings_->getString(SwooseUtilities::SettingsNames::referenceProgram);
  referenceMethod_ = settings_->getString(SwooseUtilities::SettingsNames::referenceMethod);
  referenceBasisSet_ = settings_->getString(SwooseUtilities::SettingsNames::referenceBasisSet);
  setNameOfHessianOrder();
  setNameOfStructureOptimizationOrder();
  setNameOfBondOrdersOrder();
  setNameOfAtomicChargesOrder();

  // Fill the priorities_ member with all the calculation priorities for each fragment
  evaluatePriorities();

  // Initialize the number of failed calculations counting vector
  failedCalculationsScoreForEachFragment_.resize(numberOfStructures_);
  std::fill(failedCalculationsScoreForEachFragment_.begin(), failedCalculationsScoreForEachFragment_.end(), 1);

  // Try to connect to database
  try {
    log_.output << "Connecting to database..." << Core::Log::endl;
    database_->connect();
  }
  catch (const std::exception& e) {
    throw std::runtime_error("Connecting to the specified database failed!");
  }

  // Check whether calculations collection is empty
  if (reuseDatabase_) {
    auto calculations = database_->getCollection("calculations");
    if (calculations->count("{}") == 0)
      throw std::runtime_error("Collection of calculations in given database is empty. Set \"" +
                               std::string(SwooseUtilities::SettingsNames::reuseDatabaseKey) + "\" to false.");
  }
}

DatabaseHelper::~DatabaseHelper() = default;

void DatabaseHelper::runCalculationsAndCollectResults() {
  // Resize vector of optimized structures
  data_.vectorOfOptimizedStructures.resize(numberOfStructures_);
  // Resize the vector holding the atomic charges for all subsystems
  data_.atomicChargesForEachFragment.resize(numberOfStructures_);
  // Resize the vector holding the bond order matrices for all subsystems
  data_.vectorOfBondOrderCollections.resize(numberOfStructures_);

  // Reset 'analyzed' calculations if reuseDatabase_ is true, otherwise set up initial calculations
  if (reuseDatabase_) {
    resetAnalyzedCalculationsAndDetectExistingSubsequentCalculations();
  }
  else {
    // Add structures to database
    addStructuresToDatabase();

    // Some output
    if (refineConnectivity_)
      log_.output << "Submitting structure optimizations and bond orders calculations to database..." << Core::Log::endl;
    else
      log_.output << "Submitting structure optimizations to database..." << Core::Log::endl;

    // Optimize structures and calculate bond orders for unoptimized structures
    addCalculationsToDatabaseForUnoptimizedStructures();
  }

  log_.output << "Collecting results and submitting new calculations for the optimized structures as they "
                 "become available..."
              << Core::Log::endl
              << Core::Log::endl;
  collectResultsAndSubmitSubsequentCalculations();
  log_.output << "Done." << Core::Log::endl << Core::Log::endl;
}

void DatabaseHelper::addCalculationsToDatabaseForUnoptimizedStructures() {
  auto calculations = database_->getCollection("calculations");

  // Submit structure optimization calculations
  for (int i = 0; i < numberOfStructures_; ++i) {
    DatabaseJobSubmissionHelper::submitStructureOptimization(i, calculations, structureIDs_.at(i)->string(),
                                                             priorities_.at(i), nameOfStructureOptimizationOrder_,
                                                             data_, *settings_);
  }

  // Submit bond orders calculations of original (unoptimized) fragments if required by settings
  if (refineConnectivity_) {
    for (int i = 0; i < numberOfStructures_; ++i) {
      DatabaseJobSubmissionHelper::submitBondOrdersCalculation(
          i, calculations, structureIDs_.at(i)->string(), priorities_.at(i), nameOfBondOrdersOrder_, data_, *settings_);
    }
  }
}

void DatabaseHelper::addStructuresToDatabase() {
  auto coll = database_->getCollection("structures");
  for (int i = 0; i < numberOfStructures_; ++i) {
    Database::Structure structure;
    structure.link(coll);
    int charge = data_.vectorOfChargesAndMultiplicities[i].first;
    int multiplicity = data_.vectorOfChargesAndMultiplicities[i].second;
    auto id = structure.create(*data_.vectorOfStructures[i], charge, multiplicity);
    auto label = Database::Structure::LABEL::MINIMUM_GUESS;
    structure.setLabel(label);
    structure.setComment(std::to_string(i));
    structureIDs_.push_back(std::make_unique<Database::ID>(id));
  }
}

void DatabaseHelper::collectResultsAndSubmitSubsequentCalculations() {
  // Get the collections from the database
  auto collCalc = database_->getCollection("calculations");
  auto collStruct = database_->getCollection("structures");
  auto collProp = database_->getCollection("properties");
  /*
   * The while loop will exist as long as any calculations are neither completed/analyzed nor failed.
   * The loop will also stop if enough data has been collected already.
   */
  bool calculationsMissing = true;
  FragmentDataDistributor completenessCheck(data_);
  bool enoughDataCollected = false;
  while (calculationsMissing) {
    // Verify connection to database
    while (!database_->isConnected()) {
      log_.output << "Cannot connect to database. Trying again in " << sleepTime_ << " seconds." << Core::Log::nl
                  << Core::Log::endl;
      std::this_thread::sleep_for(std::chrono::seconds(sleepTime_));
    }

    // Go through completed calculations
    int n1 = handleCompletedCalculations(collCalc, collStruct, collProp);
    // Go through failed calculations
    int n2 = handleFailedCalculations(collCalc, collStruct);

    int numberOfNewlySubmittedHessianCalculations = (n1 + n2);
    bool separateChargesCalculations = mustCalculateAtomicChargesSeparately();
    if (separateChargesCalculations)
      numberOfNewlySubmittedHessianCalculations /= 2;
    log_.debug << "Newly submitted Hessian calculations: " << numberOfNewlySubmittedHessianCalculations << Core::Log::nl;
    if (separateChargesCalculations)
      log_.debug << "Newly submitted atomic charges calculations: " << numberOfNewlySubmittedHessianCalculations
                 << Core::Log::endl;

    // Get number of incomplete calculations
    const std::string incompleteCalculationsQuery =
        "{ \"$or\" : [{ \"status\" : \"new\" }, { \"status\" : \"pending\" }] }";
    auto incompleteVec = collCalc->query<Database::Calculation>(incompleteCalculationsQuery);
    // Get number of unanalyzed calculations
    const std::string unanalyzedCalculationsQuery =
        "{ \"$or\" : [{ \"status\" : \"failed\" }, { \"status\" : \"complete\" }] }";
    auto unanalyzedVec = collCalc->query<Database::Calculation>(unanalyzedCalculationsQuery);

    // Check whether reference data is already sufficient for the parametrization
    if (earlyTerminationEnabled_ && numberOfStructures_ != 1) {
      enoughDataCollected =
          completenessCheck.referenceDataIsSufficient(refineConnectivity_, failedCalculationsScoreForEachFragment_);
      if (enoughDataCollected)
        log_.output << "Enough data has been collected." << Core::Log::endl;
    }

    // If no unfinished calculations exist, terminate the while loop, else sleep. Also terminate, if enough
    // data has been collected already and if the corresponding setting is set to true.
    if ((incompleteVec.empty() && unanalyzedVec.empty()) || enoughDataCollected)
      calculationsMissing = false;
    else {
      log_.debug << "Sleeping for " << sleepTime_ << " seconds..." << Core::Log::nl << Core::Log::endl;
      std::this_thread::sleep_for(std::chrono::seconds(sleepTime_));
    }
  }
}

int DatabaseHelper::handleCompletedCalculations(std::shared_ptr<Database::Collection> calculations,
                                                std::shared_ptr<Database::Collection> structures,
                                                std::shared_ptr<Database::Collection> properties) {
  int submitted = 0;
  const std::string queryForCompletedCalculations = "{ \"status\" : \"complete\" }";
  for (auto iter = calculations->iteratorQuery<Database::Calculation>(queryForCompletedCalculations); !iter.done(); iter++) {
    auto calc = *iter;
    calc.link(calculations);

    // Get structure
    auto structureID = calc.getStructures().at(0);
    auto structure = structures->get<Database::Structure>(structureID);
    structure.link(structures);
    int structureIndex = std::stoi(structure.getComment());

    if (structureIndex >= numberOfStructures_)
      throw std::runtime_error("At least one calculation in the database does not belong to this system.");

    // Get properties and add results to data object for different orders separately
    if (calc.getJob().order == nameOfHessianOrder_) {
      auto hessianID = structure.getProperties(nameOfHessianProperty).at(0);
      auto hessianProperty = properties->get<Database::DenseMatrixProperty>(hessianID);
      hessianProperty.link(properties);
      Utils::HessianMatrix hessian = hessianProperty.getData();
      data_.vectorOfHessians[structureIndex] = std::make_unique<Utils::HessianMatrix>(std::move(hessian));
    }
    else if (calc.getJob().order == nameOfAtomicChargesOrder_) {
      std::string nameOfChargesProperty = useGaussian_ ? nameOfCm5ChargesProperty : nameOfAtomicChargesProperty;
      auto chargesID = structure.getProperties(nameOfChargesProperty).at(0);
      auto chargesProperty = properties->get<Database::VectorProperty>(chargesID);
      chargesProperty.link(properties);
      Eigen::VectorXd chargesVector = chargesProperty.getData();
      std::vector<double> charges(chargesVector.data(), chargesVector.data() + chargesVector.size());

      // Only convert charges to CM5 if it isn't Gaussian and user requested it.
      if (!useGaussian_ && settings_->getBool(SwooseUtilities::SettingsNames::convertChargesCm5)) {
        auto convertedToCm5 =
            Utils::ChargeModel5::calculateCm5Charges(charges, *data_.vectorOfOptimizedStructures[structureIndex]);
        data_.atomicChargesForEachFragment[structureIndex] = convertedToCm5;
      }
      else {
        data_.atomicChargesForEachFragment[structureIndex] = charges;
      }
    }
    else if (calc.getJob().order == nameOfStructureOptimizationOrder_) {
      auto optStructureID = calc.getResults().structures.at(0);
      auto optStructure = structures->get<Database::Structure>(optStructureID);
      optStructure.link(structures);
      optStructure.setComment(std::to_string(structureIndex)); // forward the structure's index to the opt. structure
      data_.vectorOfOptimizedStructures[structureIndex] = std::make_unique<Utils::AtomCollection>(optStructure.getAtoms());
      if (!optimizedStructureIsValid(structureIndex))
        throw std::runtime_error(
            "At least one of the optimized structures does not match its unoptimized counterpart.");

      if (!mustCalculateAtomicChargesSeparately()) {
        std::string nameOfChargesProperty = referenceProgram_ == SwooseUtilities::OptionNames::turbomoleOption
                                                ? nameOfLoewdinChargesProperty
                                                : nameOfHirshfeldChargesProperty;
        auto chargesID = optStructure.getProperties(nameOfChargesProperty).at(0);
        auto chargesProperty = properties->get<Database::VectorProperty>(chargesID);
        chargesProperty.link(properties);
        Eigen::VectorXd chargesVector = chargesProperty.getData();
        std::vector<double> charges(chargesVector.data(), chargesVector.data() + chargesVector.size());
        // Apply CM5 correction (if requested) and store in ParametrizationData object
        if (settings_->getBool(SwooseUtilities::SettingsNames::convertChargesCm5)) {
          std::vector<double> convertedToCm5 =
              Utils::ChargeModel5::calculateCm5Charges(charges, *data_.vectorOfOptimizedStructures[structureIndex]);
          data_.atomicChargesForEachFragment[structureIndex] = convertedToCm5;
        }
        else {
          data_.atomicChargesForEachFragment[structureIndex] = charges;
        }
      }

      // Submit Hessian calculation
      bool success = DatabaseJobSubmissionHelper::submitHessianCalculation(
          structureIndex, calculations, structureID.string(), optStructureID.string(), priorities_.at(structureIndex),
          nameOfHessianOrder_, data_, *settings_, fragmentsWithHessianCalculationsInDatabase_);
      // Add to the total of the number of submitted calculations
      submitted += success ? 1 : 0;
      // Submit atomic charges calculation afterwards
      if (mustCalculateAtomicChargesSeparately()) {
        success = DatabaseJobSubmissionHelper::submitAtomicChargesCalculation(
            structureIndex, calculations, structureID.string(), optStructureID.string(), priorities_.at(structureIndex),
            nameOfAtomicChargesOrder_, data_, *settings_, fragmentsWithAtomicChargesCalculationsInDatabase_);
        // Add to the total of the number of submitted calculations
        submitted += success ? 1 : 0;
      }
    }
    else if (calc.getJob().order == nameOfBondOrdersOrder_) {
      auto bondOrdersID = structure.getProperties(nameOfBondOrdersProperty).at(0);
      auto bondOrdersProperty = properties->get<Database::SparseMatrixProperty>(bondOrdersID);
      bondOrdersProperty.link(properties);
      Eigen::SparseMatrix<double> bondOrdersMatrix = bondOrdersProperty.getData();
      Utils::BondOrderCollection bondOrders;
      bondOrders.setMatrix(std::move(bondOrdersMatrix));
      data_.vectorOfBondOrderCollections[structureIndex] =
          std::make_unique<Utils::BondOrderCollection>(std::move(bondOrders));
    }
    else {
      throw std::runtime_error("Completed calculation with illegitimate order found.");
    }

    // Modify status of calculation
    calc.setStatus(Database::Calculation::STATUS::ANALYZED);
  }
  return submitted;
}

int DatabaseHelper::handleFailedCalculations(std::shared_ptr<Database::Collection> calculations,
                                             std::shared_ptr<Database::Collection> structures) {
  int submitted = 0;
  const std::string queryForFailedCalculations = "{ \"status\" : \"failed\" }";
  for (auto iter = calculations->iteratorQuery<Database::Calculation>(queryForFailedCalculations); !iter.done(); iter++) {
    auto calc = *iter;
    calc.link(calculations);

    // Get structure
    auto structureID = calc.getStructures().at(0);
    auto structure = structures->get<Database::Structure>(structureID);
    structure.link(structures);
    int structureIndex = std::stoi(structure.getComment());

    if (structureIndex >= numberOfStructures_)
      throw std::runtime_error("At least one calculation in the database does not belong to this system.");

    // Update relevant information
    if (calc.getJob().order == nameOfHessianOrder_) {
      data_.vectorOfHessians[structureIndex] = nullptr;
      failedCalculationsScoreForEachFragment_.at(structureIndex) *= 2;
    }
    else if (calc.getJob().order == nameOfAtomicChargesOrder_) {
      data_.atomicChargesForEachFragment[structureIndex].clear(); // Should already be empty, but still for consistency
      failedCalculationsScoreForEachFragment_.at(structureIndex) *= 3;
    }
    else if (calc.getJob().order == nameOfStructureOptimizationOrder_) {
      data_.vectorOfOptimizedStructures[structureIndex] = nullptr;
      if (!mustCalculateAtomicChargesSeparately()) {
        data_.atomicChargesForEachFragment[structureIndex].clear(); // Should already be empty, but still for consistency
        failedCalculationsScoreForEachFragment_.at(structureIndex) *= 3;
      }

      // Submit Hessian calculation
      bool success = DatabaseJobSubmissionHelper::submitHessianCalculation(
          structureIndex, calculations, structureID.string(), std::string(""), priorities_.at(structureIndex),
          nameOfHessianOrder_, data_, *settings_, fragmentsWithHessianCalculationsInDatabase_);
      // Add to the total of the number of submitted calculations
      submitted += success ? 1 : 0;
      // Submit atomic charges calculation
      if (mustCalculateAtomicChargesSeparately()) {
        success = DatabaseJobSubmissionHelper::submitAtomicChargesCalculation(
            structureIndex, calculations, structureID.string(), std::string(""), priorities_.at(structureIndex),
            nameOfAtomicChargesOrder_, data_, *settings_, fragmentsWithAtomicChargesCalculationsInDatabase_);
        // Add to the total of the number of submitted calculations
        submitted += success ? 1 : 0;
      }
    }
    else if (calc.getJob().order == nameOfBondOrdersOrder_) {
      data_.vectorOfBondOrderCollections[structureIndex] = nullptr;
      failedCalculationsScoreForEachFragment_.at(structureIndex) *= 5;
    }
    else {
      throw std::runtime_error("Failed calculation with illegitimate order found.");
    }

    // Modify status of calculation
    calc.setStatus(Database::Calculation::STATUS::ANALYZED);
  }
  return submitted;
}

void DatabaseHelper::setNameOfHessianOrder() {
  nameOfHessianOrder_ = DatabaseOrderNames::nameOfScineHessianOrder;
}

void DatabaseHelper::setNameOfStructureOptimizationOrder() {
  if (referenceProgram_ == SwooseUtilities::OptionNames::orcaOption)
    nameOfStructureOptimizationOrder_ = DatabaseOrderNames::nameOfOrcaStructureOptimizationOrder;
  else if (referenceProgram_ == SwooseUtilities::OptionNames::turbomoleOption)
    nameOfStructureOptimizationOrder_ = DatabaseOrderNames::nameOfTurbomoleStructureOptimizationOrder;
  else
    nameOfStructureOptimizationOrder_ = DatabaseOrderNames::nameOfScineStructureOptimizationOrder;
}

void DatabaseHelper::setNameOfBondOrdersOrder() {
  nameOfBondOrdersOrder_ = DatabaseOrderNames::nameOfScineBondOrdersOrder;
}

void DatabaseHelper::setNameOfAtomicChargesOrder() {
  if (mustCalculateAtomicChargesSeparately())
    nameOfAtomicChargesOrder_ =
        useGaussian_ ? DatabaseOrderNames::nameOfCm5ChargesOrder : DatabaseOrderNames::nameOfScineSinglePointOrder;
}

bool DatabaseHelper::mustCalculateAtomicChargesSeparately() const {
  if (useGaussian_)
    return true;
  if (referenceProgram_ == SwooseUtilities::OptionNames::sparrowOption)
    return true;
  if (referenceProgram_ == SwooseUtilities::OptionNames::xtbOption)
    return true;
  return false;
}

void DatabaseHelper::evaluatePriorities() {
  priorities_.resize(numberOfStructures_);
  std::fill(priorities_.begin(), priorities_.end(), 10);

  if (numberOfStructures_ == 1)
    return;

  // Calculate size of largest fragment:
  int sizeOfLargestFragment = 0;
  for (int i = 0; i < numberOfStructures_; ++i) {
    auto size = data_.vectorOfStructures.at(i)->size();
    if (size > sizeOfLargestFragment)
      sizeOfLargestFragment = size;
  }

  for (int i = 0; i < numberOfStructures_; ++i) {
    // If superfluous, keep priority at 10
    if (std::find(data_.superfluousFragments.begin(), data_.superfluousFragments.end(), i) != data_.superfluousFragments.end())
      continue;

    // Case 1: Reference data generation is terminated once enough data has been collected
    if (earlyTerminationEnabled_) {
      // Gather all the priorities of the neighboring atoms
      std::vector<int> prioritiesOfNeighbors;
      for (const auto& n1 : data_.listsOfNeighbors.at(i)) {
        prioritiesOfNeighbors.push_back(priorities_.at(n1));
        for (const auto& n2 : data_.listsOfNeighbors.at(n1)) {
          if (n2 == i)
            continue;
          prioritiesOfNeighbors.push_back(priorities_.at(n2));
        }
      }

      /*
       * If no priority of 1 or 2 exists yet in the environment, give it to this atom, otherwise increase priority
       * iteratively and proceed with the same logic.
       * The exact priority is determined based on the fragment size. Fragments that have at least 70%
       * of the maximum number of atoms per fragment, get a high priority.
       */
      for (int p = 1; p < 10; p += 2) {
        int q = p + 1;
        if (std::find(prioritiesOfNeighbors.begin(), prioritiesOfNeighbors.end(), p) != prioritiesOfNeighbors.end())
          continue;
        if (std::find(prioritiesOfNeighbors.begin(), prioritiesOfNeighbors.end(), q) != prioritiesOfNeighbors.end())
          continue;

        // Set this priority to either p or q based on the fragment size
        if (data_.vectorOfStructures.at(i)->size() > 0.7 * sizeOfLargestFragment)
          priorities_.at(i) = p;
        else
          priorities_.at(i) = q;
        break;
      }
    }
    else { // Case 2: Reference data generation is not terminated once enough data has been collected
      // Priority just based on size of system, project size of fragment to interval 1 to 10.
      double p = (data_.vectorOfStructures.at(i)->size() * 9.0 / sizeOfLargestFragment) + 1.0;
      priorities_.at(i) = 11 - static_cast<int>(std::round(p));
    }
  }
}

void DatabaseHelper::resetAnalyzedCalculationsAndDetectExistingSubsequentCalculations() {
  auto calculations = database_->getCollection("calculations");
  auto structures = database_->getCollection("structures");
  const std::string emptyQuery = "{}";
  for (auto iter = calculations->iteratorQuery<Database::Calculation>(emptyQuery); !iter.done(); iter++) {
    auto calc = *iter;
    calc.link(calculations);

    // Detect existing Hessian and atomic charges calculations
    bool isHessianCalc = calc.getJob().order == nameOfHessianOrder_;
    bool isAtomicChargesCalc = calc.getJob().order == nameOfAtomicChargesOrder_;
    if (isHessianCalc || isAtomicChargesCalc) {
      // Get fragment index
      auto structureID = calc.getStructures().at(0);
      auto structure = structures->get<Database::Structure>(structureID);
      structure.link(structures);
      int structureIndex = std::stoi(structure.getComment());
      if (isHessianCalc)
        fragmentsWithHessianCalculationsInDatabase_.insert(structureIndex);
      else
        fragmentsWithAtomicChargesCalculationsInDatabase_.insert(structureIndex);
    }

    // Reset calculation status if it is 'analyzed'
    if (calc.getStatus() == Database::Calculation::STATUS::ANALYZED) {
      bool hasFailed = analyzedCalculationsHasFailed(calc);
      if (hasFailed)
        calc.setStatus(Database::Calculation::STATUS::FAILED);
      else
        calc.setStatus(Database::Calculation::STATUS::COMPLETE);
    }
  }
}

bool DatabaseHelper::analyzedCalculationsHasFailed(const Database::Calculation& calculation) const {
  if (calculation.getStatus() != Database::Calculation::STATUS::ANALYZED)
    throw std::runtime_error(
        "The calculation that was passed to function 'analyzedCalculationsHasFailed()' has not been analyzed yet.");

  bool isStructureOptimization = calculation.getJob().order == nameOfStructureOptimizationOrder_;
  if (isStructureOptimization) {
    return calculation.getResults().structures.empty();
  }
  return calculation.getResults().properties.empty();
}

bool DatabaseHelper::optimizedStructureIsValid(int structureIndex) const {
  if (data_.vectorOfOptimizedStructures.at(structureIndex) == nullptr || data_.vectorOfStructures.at(structureIndex) == nullptr)
    return false;
  if (data_.vectorOfOptimizedStructures.at(structureIndex)->size() != data_.vectorOfStructures.at(structureIndex)->size())
    return false;
  for (int i = 0; i < data_.vectorOfOptimizedStructures.at(structureIndex)->size(); ++i) {
    Utils::ElementType e1 = data_.vectorOfOptimizedStructures.at(structureIndex)->getElement(i);
    Utils::ElementType e2 = data_.vectorOfStructures.at(structureIndex)->getElement(i);
    if (e1 != e2)
      return false;
  }
  return true;
}

void DatabaseHelper::dropDatabase() {
  if (database_->isConnected())
    database_->wipe();
}

} // namespace MMParametrization
} // namespace Scine
