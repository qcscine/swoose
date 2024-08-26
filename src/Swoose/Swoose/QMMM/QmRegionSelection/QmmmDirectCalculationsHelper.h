/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMMMDIRECTCALCULATIONSHELPER_H
#define SWOOSE_QMMM_QMMMDIRECTCALCULATIONSHELPER_H

#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace Utils {
class Settings;
class AtomCollection;
class BondOrderCollection;
} // namespace Utils

namespace Qmmm {
struct QmmmModel;
struct QmmmData;
class QmmmCalculator;
using ForcesCollection = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

/**
 * @class QmmmDirectCalculationsHelper QmmmDirectCalculationsHelper.h
 * @brief This class handles QM/MM calculations in the direct mode.
 */
class QmmmDirectCalculationsHelper {
 public:
  /**
   * @brief Constructor.
   */
  QmmmDirectCalculationsHelper(std::shared_ptr<QmmmCalculator> qmmmCalculator, const Utils::Settings& settings, Core::Log& log,
                               const Utils::AtomCollection& structure, const std::vector<QmmmModel>& qmmmModelCandidates,
                               const std::vector<QmmmModel>& qmmmReferenceModels, const QmmmData& qmmmData);

  /**
   * @brief Calculates the forces for the QM/MM candidate and reference models.
   * @return The forces obtained with the QM/MM candidate and reference models.
   */
  std::vector<ForcesCollection> calculateForces();

 private:
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
  // The data calculated for the QM/MM candidate and reference models.
  const QmmmData& qmmmData_;
  std::shared_ptr<Qmmm::QmmmCalculator> qmmmCalculator_;
};

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMMMDIRECTCALCULATIONSHELPER_H
