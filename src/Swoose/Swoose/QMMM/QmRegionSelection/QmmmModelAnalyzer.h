/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMMMMODELANALYZER_H
#define SWOOSE_QMMM_QMMMMODELANALYZER_H

#include <Utils/Constants.h>
#include <Eigen/Dense>
#include <vector>

namespace Scine {

namespace Core {
class Log;
} // namespace Core

namespace Utils {
class Settings;
class AtomCollection;
} // namespace Utils

namespace Qmmm {
struct QmmmData;
struct QmmmModel;

/**
 * @class QmmmModelAnalyzer QmmmModelAnalyzer.h
 * @brief TODO
 */
class QmmmModelAnalyzer {
 public:
  /**
   * @brief Constructor. It also already performs the analysis based on the given reference data.
   * @param settings The settings of the QmRegionSelector.
   * @param log The logger.
   * @param data The reference data.
   * @param structure The molecular structure of the system.
   * @param qmmmModelCandidates A vector containing the QM/MM model candidates.
   */
  QmmmModelAnalyzer(const Utils::Settings& settings, Core::Log& log, const QmmmData& data,
                    const Utils::AtomCollection& structure, const std::vector<QmmmModel>& qmmmModelCandidates);

  /**
   * @brief Returns the index of the candidate model that is the optimal choice based on the reference data.
   * @return The index of the candidate model that is the optimal choice based on the reference data.
   */
  int getIndexOfOptimalModel() const;

 private:
  // Function that performs the analysis of the reference data
  void analyzeData();
  // Function returning a vector of atom indices of the atoms that are close to the center atom
  std::vector<int> getAtomIndicesCloseToCenterAtom();
  std::vector<Eigen::RowVector3d> calculateReferenceForces(const std::vector<int>& relevantAtoms);
  // Calculates the mean error of the forces for a given candidate model (with index 'modelIndex')
  double calculateMeanErrorForCandidateModel(int modelIndex, const std::vector<int>& relevantAtoms,
                                             const std::vector<Eigen::RowVector3d>& referenceForces);
  // The settings.
  const Utils::Settings& settings_;
  // Logger.
  Core::Log& log_;
  // The QM/MM reference data
  const QmmmData& data_;
  // Molecular structure of the whole system
  const Utils::AtomCollection& structure_;
  // QM/MM model candidates
  const std::vector<QmmmModel>& candidates_;
  // Index of the optimal candidate model, default of -1 signals that no decision has been made yet
  int indexOfOptimalModel_ = -1;
  static constexpr double distanceThresholdForAnalysis_ = 4.0 * Utils::Constants::bohr_per_angstrom; // 4 Angstrom
};

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMMMMODELANALYZER_H
