/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_DATABASEHELPER_H
#define MMPARAMETRIZATION_DATABASEHELPER_H

#include <memory>
#include <unordered_set>
#include <vector>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace Utils {
class Settings;
} // namespace Utils

namespace Database {
class Manager;
class ID;
class Collection;
class Calculation;
} // namespace Database

namespace MMParametrization {
struct ParametrizationData;

class DatabaseHelper {
 public:
  /**
   * @brief Constructor.
   */
  DatabaseHelper(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings, Core::Log& log);
  /**
   * @brief Destructor.
   */
  ~DatabaseHelper();
  /**
   * @brief Deletes the database.
   */
  void dropDatabase();
  /**
   * @brief Runs structure optimizations and other reference calculations (Hessians, charges, bond orders)
   *        and saves the results in the corresponding fields in the ParametrizationData object.
   */
  void runCalculationsAndCollectResults();

 private:
  /*
   * @brief Resets all 'analyzed' calculations in the database to either 'complete' or 'failed' depending on
   *        the present results. Furthermore, this function detects and stores which of the fragments already have
   *        their Hessian calculations (and if needed atomic charges) in the database.
   */
  void resetAnalyzedCalculationsAndDetectExistingSubsequentCalculations();
  /*
   * @brief Returns whether a given analyzed calculation has failed.
   * @throws std::runtime_error If given calculations does not have the status 'analyzed'.
   */
  bool analyzedCalculationsHasFailed(const Database::Calculation& calculation) const;
  /*
   * @brief This function adds the structure optimization calculations and the bond orders calculations
   *        to the corresponding collection in the database.
   */
  void addCalculationsToDatabaseForUnoptimizedStructures();
  /*
   * @brief This function collects the results of the reference calculations once they are finished.
   *        Furthermore, it submits Hessian and atomic charges calculations if a structure optimizations is done.
   */
  void collectResultsAndSubmitSubsequentCalculations();
  /*
   * @brief Iterates over all completed calculations in the database and collects the results.
   *        Furthermore, it submits Hessian and atomic charges calculations if a structure optimizations is done.
   *
   * @return The number of newly submitted calculations.
   */
  int handleCompletedCalculations(std::shared_ptr<Database::Collection> calculations,
                                  std::shared_ptr<Database::Collection> structures,
                                  std::shared_ptr<Database::Collection> properties);
  /*
   * @brief Iterates over all failed calculations in the database and handles those cases.
   *        Furthermore, it submits Hessian and atomic charges calculations if a structure optimizations is done.
   *
   * @return The number of newly submitted calculations.
   */
  int handleFailedCalculations(std::shared_ptr<Database::Collection> calculations,
                               std::shared_ptr<Database::Collection> structures);

  // Adds the structures of the subsystems to the corresponding collection in the database.
  void addStructuresToDatabase();
  /*
   * @brief Internally sets the correct name of the Hessian order.
   */
  void setNameOfHessianOrder();
  /*
   * @brief Internally sets the correct name of the structure optimization order.
   */
  void setNameOfStructureOptimizationOrder();
  /*
   * @brief Internally sets the correct name of the bond orders order.
   */
  void setNameOfBondOrdersOrder();
  /*
   * @brief Internally sets the correct name of the atomic charges order.
   */
  void setNameOfAtomicChargesOrder();
  /*
   * @brief Decides whether atomic charges must be calculated in a separate calculation.
   */
  bool mustCalculateAtomicChargesSeparately() const;
  /*
   * @brief Fills the vector of priorities for all the calculations (one priority per fragment).
   */
  void evaluatePriorities();
  /*
   * @brief Returns whether the optimized structure with the given
   *        structure index is valid compared to the unoptimized one.
   */
  bool optimizedStructureIsValid(int structureIndex) const;
  /*
   * Vector of database IDs corresponding to the structures in the database.
   */
  std::vector<std::unique_ptr<Database::ID>> structureIDs_;
  // Number of subsystems.
  int numberOfStructures_;
  // Sleep time in between database operations in seconds
  int sleepTime_;
  // Whether to calculate bond orders to refine the initial connectivity.
  bool refineConnectivity_;
  // Whether to use Gaussian as a program for calculating the CM5 charges.
  bool useGaussian_;
  // Whether to terminate the reference data generation early when enough data is collected
  bool earlyTerminationEnabled_;
  /*
   * Tracks the failed calculations for each fragment in a vector in the following way:
   * - multiple of 2 if Hessian failed
   * - multiple of 3 if atomic charges failed
   * - multiple of 5 if bond orders failed
   *
   * This allows to reconstruct which calculations already failed later, to decide whether the parametrization
   * cannot be completed anymore. Whether or not structure optimization failed is not important to track because
   * if a structure optimization fails, then the Hessian also fails immediately.
   */
  std::vector<int> failedCalculationsScoreForEachFragment_;
  // Whether an existing database is exploited
  bool reuseDatabase_;
  // In the case of reusing an existing database, this keeps track of existing Hessian calculations.
  std::unordered_set<int> fragmentsWithHessianCalculationsInDatabase_;
  // In the case of reusing an existing database, this keeps track of existing atomic charges calculations.
  std::unordered_set<int> fragmentsWithAtomicChargesCalculationsInDatabase_;
  // Reference Program
  std::string referenceProgram_;
  // Reference Method
  std::string referenceMethod_;
  // Reference Basis Set
  std::string referenceBasisSet_;
  // Name of the structure optimization order
  std::string nameOfStructureOptimizationOrder_;
  // Name of the Hessian calculation order
  std::string nameOfHessianOrder_;
  // Name of the bond orders order
  std::string nameOfBondOrdersOrder_;
  // Name of the atomic charges order (empty if no such calculations are performed)
  std::string nameOfAtomicChargesOrder_;
  // The vector of priorities of each fragment
  std::vector<int> priorities_;
  // The data used within all MM parametrization classes
  ParametrizationData& data_;
  // The settings
  std::shared_ptr<Utils::Settings> settings_;
  // The database
  std::unique_ptr<Database::Manager> database_;
  // The logger.
  Core::Log& log_;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_DATABASEHELPER_H
