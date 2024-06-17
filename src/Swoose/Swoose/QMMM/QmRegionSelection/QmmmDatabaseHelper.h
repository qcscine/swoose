/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMMMDATABASEHELPER_H
#define SWOOSE_QMMM_QMMMDATABASEHELPER_H

#include <Database/Manager.h>
#include <Utils/IO/Yaml.h>
#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Scine {

namespace Core {
class Log;
} // namespace Core

namespace Utils {
class Settings;
class AtomCollection;
class BondOrderCollection;
} // namespace Utils

namespace Database {
class ID;
class Collection;
} // namespace Database

namespace Qmmm {
struct QmmmModel;
struct QmmmData;
using ForcesCollection = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

/**
 * @class QmmmDatabaseHelper QmmmDatabaseHelper.h
 * @brief This class handles the submission and collection of results for QM/MM calculations via the database.
 */
class QmmmDatabaseHelper {
 public:
  /**
   * @brief Constructor.
   */
  QmmmDatabaseHelper(const Utils::Settings& settings, Core::Log& log, const Utils::AtomCollection& structure,
                     const Utils::BondOrderCollection& bondOrders, const std::vector<QmmmModel>& qmmmModelCandidates,
                     const std::vector<QmmmModel>& qmmmReferenceModels, const QmmmData& qmmmData,
                     const Utils::Settings& calculatorSettings);

  /**
   * @brief Calculates the forces for the QM/MM candidate and reference models.
   * @return The forces obtained with the QM/MM candidate and reference models.
   */
  std::vector<ForcesCollection> calculateForces();

 private:
  // Returns the content of the parameters file as a string
  std::string getParametersAsString() const;
  // Submit a QM/MM reference calculation
  bool submitCalculation(const QmmmModel& qmmmModel, const Database::ID& structureID, int qmmmModelIndex,
                         std::shared_ptr<Database::Collection> calculations, double symmetryScore = 0.0,
                         double maxAllowedSymmetryScore = 0.0);
  // Collect the results and return the results
  std::vector<ForcesCollection> collectCalculationResults(std::shared_ptr<Database::Collection> calculations,
                                                          std::shared_ptr<Database::Collection> structures,
                                                          std::shared_ptr<Database::Collection> properties);
  // Logger.
  Core::Log& log_;
  // The settings.
  const Utils::Settings& settings_;
  // The database
  std::unique_ptr<Database::Manager> database_;
  // QM/MM model candidates
  const std::vector<QmmmModel>& qmmmModelCandidates_;
  // QM/MM reference models
  const std::vector<QmmmModel>& qmmmReferenceModels_;
  // Molecular structure of the whole system
  const Utils::AtomCollection& structure_;
  // Sleep time in between database operations in seconds
  int sleepTime_;
  // The bond orders of the system
  const Utils::BondOrderCollection& bondOrders_;
  // The data calculated for the QM/MM candidate and reference models.
  const QmmmData& qmmmData_;
  // The settings of the QM/MM calculator
  const Utils::Settings& calculatorSettings_;
};

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMMMDATABASEHELPER_H
