/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DatabaseJobSubmissionHelper.h"
#include "../CalculationManager.h"
#include "../MMParametrizationSettings.h"
#include "../ParametrizationData.h"
#include "BasicJobSubmissionHelper.h"
#include "DatabaseOrderNames.h"
#include <Database/Collection.h>
#include <Database/Manager.h>
#include <Database/Objects/Calculation.h>
#include <Database/Objects/DenseMatrixProperty.h>
#include <Database/Objects/Structure.h>

namespace Scine {
namespace MMParametrization {
namespace DatabaseJobSubmissionHelper {

void submitStructureOptimization(int fragmentIndex, std::shared_ptr<Database::Collection> calcsColl,
                                 std::string structureIDString, int priority, const std::string& orderName,
                                 const ParametrizationData& data, const Utils::Settings& settings) {
  Database::ID structureID(structureIDString);
  Database::Calculation calc;
  calc.link(calcsColl);

  // Generate the job specifications
  Database::Calculation::Job job(orderName);
  job.cores = settings.getInt(Utils::SettingsNames::externalProgramNProcs);
  job.memory = 1.0 + (job.cores / 2) + (data.vectorOfStructures.at(fragmentIndex)->size() / 10.0);
  job.disk = 4.0 + (data.vectorOfStructures.at(fragmentIndex)->size() / 10.0);

  auto method = settings.getString(SwooseUtilities::SettingsNames::referenceMethod);
  auto referenceProgram = settings.getString(SwooseUtilities::SettingsNames::referenceProgram);
  Database::Model model(BasicJobSubmissionHelper::determineMethodFamily(method, referenceProgram), method,
                        BasicJobSubmissionHelper::determineBasisSet(settings));
  model.program = referenceProgram;
  calc.create(model, job, {structureID});

  // Add Cartesian constraints
  const auto constrainedAtoms = data.constrainedAtoms.at(fragmentIndex);

  // Set constraints // TODO: Unify the settings names eventually.
  if (!constrainedAtoms.empty()) {
    if (orderName == DatabaseOrderNames::nameOfScineStructureOptimizationOrder) {
      calc.setSetting("geoopt_coordinate_system", "cartesian");
      calc.setSetting("geoopt_constrained_atoms", constrainedAtoms);
    }
    else {
      calc.setSetting("cartesian_constraints", constrainedAtoms);
    }
  }

  calc.setSetting("convergence_max_iterations", 3000);           // TODO: should this always be so large?
  calc.setSetting(Utils::SettingsNames::maxScfIterations, 1000); // TODO: should this always be so large?
  applyScfSafetySettings(calc, settings);
  // Important in particular for SCINE optimizations, but also helpful for the others:
  calc.setSetting(Utils::SettingsNames::selfConsistenceCriterion, 1e-7);

  if (orderName == DatabaseOrderNames::nameOfScineStructureOptimizationOrder) {
    calc.setSetting("optimizer", "bfgs");
    calc.setSetting("bfgs_use_trust_radius", true);
  }

  if (!settings.getBool(SwooseUtilities::SettingsNames::useGaussianOptionKey)) {
    if (referenceProgram == SwooseUtilities::OptionNames::turbomoleOption)
      calc.setSetting("calculate_loewdin_charges", true);
    else if (referenceProgram == SwooseUtilities::OptionNames::orcaOption)
      calc.setSetting("calculate_hirshfeld_charges", true);
  }

  // Set priority of calculation
  calc.setPriority(priority);

  // Once the status is set to "new", puffin can start fetching it.
  // Set the status to failed directly if the fragment is superfluous.
  if (std::find(data.superfluousFragments.begin(), data.superfluousFragments.end(), fragmentIndex) ==
      data.superfluousFragments.end())
    calc.setStatus(Database::Calculation::STATUS::NEW);
  else
    calc.setStatus(Database::Calculation::STATUS::FAILED);
}

void submitBondOrdersCalculation(int fragmentIndex, std::shared_ptr<Database::Collection> calcsColl,
                                 std::string structureIDString, int priority, const std::string& orderName,
                                 const ParametrizationData& data, const Utils::Settings& settings) {
  Database::ID structureID(structureIDString);
  Database::Calculation calc;
  calc.link(calcsColl);
  // Generate the job specifications
  Database::Calculation::Job job(orderName);
  job.memory = 1.0 + (data.vectorOfStructures.at(fragmentIndex)->size() / 20.0);
  job.cores = 1;
  job.disk = 4.0 + (data.vectorOfStructures.at(fragmentIndex)->size() / 10.0);
  auto method = settings.getString(SwooseUtilities::SettingsNames::referenceMethod);
  auto referenceProgram = settings.getString(SwooseUtilities::SettingsNames::referenceProgram);
  Database::Model model(BasicJobSubmissionHelper::determineMethodFamily(method, referenceProgram), method,
                        BasicJobSubmissionHelper::determineBasisSet(settings));
  model.program = referenceProgram;
  calc.create(model, job, {structureID});

  calc.setSetting(Utils::SettingsNames::maxScfIterations, 1000); // TODO: should this always be so large?
  applyScfSafetySettings(calc, settings);

  // Set priority of calculation
  calc.setPriority(priority);

  // Once the status is set to "new", puffin can start fetching it.
  // Set the status to failed directly if the fragment is superfluous.
  if (std::find(data.superfluousFragments.begin(), data.superfluousFragments.end(), fragmentIndex) ==
      data.superfluousFragments.end())
    calc.setStatus(Database::Calculation::STATUS::NEW);
  else
    calc.setStatus(Database::Calculation::STATUS::FAILED);
}

bool submitHessianCalculation(int fragmentIndex, std::shared_ptr<Database::Collection> calcsColl,
                              std::string unoptimizedStructureIDString, std::string optimizedStructureIDString, int priority,
                              const std::string& orderName, const ParametrizationData& data, const Utils::Settings& settings,
                              const std::unordered_set<int>& fragmentsWithHessianCalculationsInDatabase) {
  // Only submit Hessian calculation if it isn't already in the database (during reusing database run)
  if (fragmentsWithHessianCalculationsInDatabase.find(fragmentIndex) != fragmentsWithHessianCalculationsInDatabase.end())
    return false;
  Database::Calculation calc;
  calc.link(calcsColl);

  // Generate the job specifications
  Database::Calculation::Job job(orderName);
  job.cores = settings.getInt(Utils::SettingsNames::externalProgramNProcs);
  job.memory = 1.0 + (job.cores / 2) + (data.vectorOfStructures.at(fragmentIndex)->size() / 10.0);
  job.disk = 4.0 + (data.vectorOfStructures.at(fragmentIndex)->size() / 10.0);

  auto method = settings.getString(SwooseUtilities::SettingsNames::referenceMethod);
  auto referenceProgram = settings.getString(SwooseUtilities::SettingsNames::referenceProgram);
  Database::Model model(BasicJobSubmissionHelper::determineMethodFamily(method, referenceProgram), method,
                        BasicJobSubmissionHelper::determineBasisSet(settings));
  model.program = referenceProgram;

  if (!optimizedStructureIDString.empty()) {
    Database::ID optimizedStructureID(optimizedStructureIDString);
    calc.create(model, job, {optimizedStructureID});
    // Pass on the internal settings for the calculation to the database
    calc.setSetting(Utils::SettingsNames::maxScfIterations, 1000); // TODO: should this always be so large?
    applyScfSafetySettings(calc, settings);
    calc.setSetting(Utils::SettingsNames::selfConsistenceCriterion, 1e-7);

    // Set priority of calculation
    calc.setPriority(priority);

    // Once the status is set to "new", puffin can start fetching it.
    calc.setStatus(Database::Calculation::STATUS::NEW);
    return true;
  }
  // Create a calculation with the unoptimized structure to retain the fragment index
  Database::ID unoptimizedStructureID(unoptimizedStructureIDString);
  calc.create(model, job, {unoptimizedStructureID});

  // Set priority of calculation
  calc.setPriority(priority);

  // If the optimization failed, this is also an automatic fail.
  calc.setStatus(Database::Calculation::STATUS::FAILED);
  return false;
}

bool submitAtomicChargesCalculation(int fragmentIndex, std::shared_ptr<Database::Collection> calcsColl,
                                    std::string unoptimizedStructureIDString, std::string optimizedStructureIDString,
                                    int priority, const std::string& orderName, const ParametrizationData& data,
                                    const Utils::Settings& settings,
                                    const std::unordered_set<int>& fragmentsWithAtomicChargesCalculationsInDatabase) {
  // Only submit atomic charges calculation if it isn't already in the database (during reusing database run)
  if (fragmentsWithAtomicChargesCalculationsInDatabase.find(fragmentIndex) !=
      fragmentsWithAtomicChargesCalculationsInDatabase.end())
    return false;
  Database::Calculation calc;
  calc.link(calcsColl);
  // Generate the job specifications
  Database::Calculation::Job job(orderName);
  job.memory = 1.0 + (data.vectorOfStructures.at(fragmentIndex)->size() / 20.0);
  job.cores = 1;
  job.disk = 4.0 + (data.vectorOfStructures.at(fragmentIndex)->size() / 10.0);

  Database::Model model("", "", "");
  if (settings.getBool(SwooseUtilities::SettingsNames::useGaussianOptionKey)) {
    model.methodFamily = "dft";
    model.method = settings.getString(SwooseUtilities::SettingsNames::gaussianMethod);
    model.basisSet = settings.getString(SwooseUtilities::SettingsNames::gaussianBasisSet);
    model.program = "gaussian";
  }
  else {
    auto referenceProgram = settings.getString(SwooseUtilities::SettingsNames::referenceProgram);
    model.method = settings.getString(SwooseUtilities::SettingsNames::referenceMethod);
    model.methodFamily = BasicJobSubmissionHelper::determineMethodFamily(model.method, referenceProgram);
    model.basisSet = BasicJobSubmissionHelper::determineBasisSet(settings);
    model.program = referenceProgram;
  }

  if (!optimizedStructureIDString.empty()) {
    Database::ID optimizedStructureID(optimizedStructureIDString);
    calc.create(model, job, {optimizedStructureID});

    // Set priority of calculation
    calc.setPriority(priority);

    // Once the status is set to "new", puffin can start fetching it.
    calc.setStatus(Database::Calculation::STATUS::NEW);
    return true;
  }

  // Create a calculation with the unoptimized structure to retain the fragment index
  Database::ID unoptimizedStructureID(unoptimizedStructureIDString);
  calc.create(model, job, {unoptimizedStructureID});

  // Set priority of calculation
  calc.setPriority(priority);

  // If the optimization failed, this is also an automatic fail.
  calc.setStatus(Database::Calculation::STATUS::FAILED);
  return false;
}

void applyScfSafetySettings(Database::Calculation& calculation, const Utils::Settings& settings) {
  if (!settings.getBool(SwooseUtilities::SettingsNames::increaseScfSafetyKey))
    return;

  std::string referenceProgram = settings.getString(SwooseUtilities::SettingsNames::referenceProgram);
  if (referenceProgram == SwooseUtilities::OptionNames::turbomoleOption) {
    calculation.setSetting("scf_orbitalshift", 0.5);
    calculation.setSetting("scf_damping", true);
  }
  else if (referenceProgram == SwooseUtilities::OptionNames::orcaOption) {
    calculation.setSetting("scf_damping", true);
  }
}

} // namespace DatabaseJobSubmissionHelper
} // namespace MMParametrization
} // namespace Scine
