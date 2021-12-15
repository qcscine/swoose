/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_SYMMETRYSCORES_H
#define SWOOSE_QMMM_SYMMETRYSCORES_H

#include <vector>

namespace Scine {

namespace Utils {
class Settings;
} // namespace Utils

namespace Qmmm {
struct QmmmModel;
struct QmmmData;

/**
 * @brief Calculates the symmetry score for one QM/MM model.
 * @param model The model.
 * @param distances Vector containing the distances of all atoms to the center atom of the QM region.
 * @param numberOfValuesForMean The number of values X that are taken to calculate the mean distance of the
 *                              X most distant QM atoms (from the center) and the X least distant non-QM atoms.
 * @return The symmetry score. If there are less than 'numberOfValuesForMean' QM or non-QM atoms, a symmetry score
 *         of 0.0 is returned.
 */
double calculateSymmetryScoreForOneModel(const QmmmModel& model, const std::vector<double>& distances, int numberOfValuesForMean);

/**
 * @brief Calculates the maximum allowed symmetry score based on the symmetry scores of all the candidate models,
 *        which can be found in the QmmmData object.
 * @brief The reference data of all QM/MM structure separations. Here, we need the calculated symmetry scores.
 * @brief The settings of the QM region selection.
 * @return The maximum allowed symmetry score based on the symmetry scores of all the candidate models.
 */
double calculateMaxAllowedSymmetryScore(const QmmmData& qmmmData, const Utils::Settings& settings);

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMMMREFERENCEDATAMANAGER_H
