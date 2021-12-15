/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMMMREFERENCEDATAMANAGER_H
#define SWOOSE_QMMM_QMMMREFERENCEDATAMANAGER_H

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

namespace Qmmm {
struct QmmmModel;
using ForcesCollection = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

/**
 * @struct QmmmReferenceData QmmmReferenceDataManager.h
 * @brief Container for the reference data for all of the QM/MM models.
 */
struct QmmmData {
  /// @brief Container for the atomic forces calculated with each QM/MM model.
  std::vector<ForcesCollection> forces;
  /// @brief Container for the symmetry scores for each QM/MM model.
  std::vector<double> symmetryScores;
  /// @brief Container for the number of link atoms for each QM/MM model.
  std::vector<int> linkAtomNumbers;
  /// @brief Number of reference calculations. The last nRef data points in the objects above correspond to these.
  int nRef = 0;
};

/**
 * @class QmmmReferenceDataManager QmmmReferenceDataManager.h
 * @brief Manages the calculation of all data contained in the QmmmData struct.
 */
class QmmmReferenceDataManager {
 public:
  /**
   * @brief Constructor.
   */
  QmmmReferenceDataManager(const Utils::Settings& settings, Core::Log& log, const Utils::AtomCollection& structure,
                           const Utils::BondOrderCollection& bondOrders, const std::vector<QmmmModel>& qmmmModelCandidates,
                           const std::vector<QmmmModel>& qmmmReferenceModels);

  /**
   * @brief Calculates the properties in the QmmmData struct for the QM/MM candidate models and reference models.
   * @return The data that was calculated with the QM/MM candidate models and reference models.
   */
  QmmmData calculateData();

 private:
  // Calculates the symmetry scores and stores them in the QM/MM data
  void calculateSymmetryScores(QmmmData& data) const;
  // Calculates the number of link atoms for each QM/MM model and stores these values in the QM/MM data
  void calculateLinkAtomNumbers(QmmmData& data) const;
  // Perform reference calculations and store results in QmmmData object
  void handleReferenceCalculations(QmmmData& data);
  // The settings.
  const Utils::Settings& settings_;
  // Logger.
  Core::Log& log_;
  // Molecular structure of the whole system
  const Utils::AtomCollection& structure_;
  // QM/MM model candidates
  const std::vector<QmmmModel>& qmmmModelCandidates_;
  // QM/MM reference models
  const std::vector<QmmmModel>& qmmmReferenceModels_;
  // The bond orders of the system
  const Utils::BondOrderCollection& bondOrders_;
};

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMMMREFERENCEDATAMANAGER_H
