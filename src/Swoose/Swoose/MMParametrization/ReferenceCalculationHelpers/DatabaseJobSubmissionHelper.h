/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_DATABASEJOBSUBMISSIONHELPER_H
#define MMPARAMETRIZATION_DATABASEJOBSUBMISSIONHELPER_H

#include <memory>
#include <unordered_set>
#include <vector>

namespace Scine {

namespace Utils {
class Settings;
} // namespace Utils

namespace Database {
class ID;
class Collection;
class Calculation;
} // namespace Database

namespace MMParametrization {
struct ParametrizationData;

namespace DatabaseJobSubmissionHelper {

/**
 * @brief This function adds one structure optimization to the database (given collection)
 *        for the fragment with the given index.
 * @param fragmentIndex Index of the fragment.
 * @param calcsColl Calculations collection in the database.
 * @param structureIDString String representation of database ID for the fragment structure.
 * @param priority Priority of the calculation. 1 is the highest, 10 the lowest.
 * @param orderName Name of the structure optimization order.
 * @param data The ParametrizationData object.
 * @param settings The settings of the parametrization.
 */
void submitStructureOptimization(int fragmentIndex, std::shared_ptr<Database::Collection> calcsColl,
                                 std::string structureIDString, int priority, const std::string& orderName,
                                 const ParametrizationData& data, const Utils::Settings& settings);

/**
 * @brief This function adds one bond orders calculation to the database (given collection)
 *        for the fragment with the given index.
 * @param fragmentIndex Index of the fragment.
 * @param calcsColl Calculations collection in the database.
 * @param structureIDString String representation of database ID for the fragment structure.
 * @param priority Priority of the calculation. 1 is the highest, 10 the lowest.
 * @param orderName Name of the bond orders order.
 * @param data The ParametrizationData object.
 * @param settings The settings of the parametrization.
 */
void submitBondOrdersCalculation(int fragmentIndex, std::shared_ptr<Database::Collection> calcsColl,
                                 std::string structureIDString, int priority, const std::string& orderName,
                                 const ParametrizationData& data, const Utils::Settings& settings);

/**
 * @brief This function adds one Hessian calculation to the database (given collection)
 *        for the fragment with the given index.
 * @param fragmentIndex Index of the fragment.
 * @param calcsColl Calculations collection in the database.
 * @param unoptimizedStructureIDString String representation of database ID for the fragment structure.
 *                                     Empty if not available.
 * @param optimizedStructureIDString String representation of database ID for the fragment structure.
 *                                   Empty if not available.
 * @param priority Priority of the calculation. 1 is the highest, 10 the lowest.
 * @param orderName Name of the Hessian order.
 * @param data The ParametrizationData object.
 * @param settings The settings of the parametrization.
 * @param fragmentsWithHessianCalculationsInDatabase Set of fragment indices that already have Hessian calculations
 *        in the database.
 * @return Whether a new (not directly failed) calculation was put into the database.
 */
bool submitHessianCalculation(int fragmentIndex, std::shared_ptr<Database::Collection> calcsColl,
                              std::string unoptimizedStructureIDString, std::string optimizedStructureIDString, int priority,
                              const std::string& orderName, const ParametrizationData& data, const Utils::Settings& settings,
                              const std::unordered_set<int>& fragmentsWithHessianCalculationsInDatabase);
/**
 * @brief This function adds one atomic charges calculation to the database (given collection)
 *        for the fragment with the given index.
 * @param fragmentIndex Index of the fragment.
 * @param calcsColl Calculations collection in the database.
 * @param unoptimizedStructureIDString String representation of database ID for the fragment structure.
 *                                     Empty if not available.
 * @param optimizedStructureIDString String representation of database ID for the fragment structure.
 *                                   Empty if not available.
 * @param priority Priority of the calculation. 1 is the highest, 10 the lowest.
 * @param orderName Name of the atomic charges order.
 * @param data The ParametrizationData object.
 * @param settings The settings of the parametrization.
 * @param fragmentsWithAtomicChargesCalculationsInDatabase Set of fragment indices that already have atomic charges
 *        calculations in the database.
 * @return Whether a new (not directly failed) calculation was put into the database.
 */
bool submitAtomicChargesCalculation(int fragmentIndex, std::shared_ptr<Database::Collection> calcsColl,
                                    std::string unoptimizedStructureIDString, std::string optimizedStructureIDString,
                                    int priority, const std::string& orderName, const ParametrizationData& data,
                                    const Utils::Settings& settings,
                                    const std::unordered_set<int>& fragmentsWithAtomicChargesCalculationsInDatabase);

/**
 * @brief If possible and requested, adds setting for more SCF convergence safety to the given calculation.
 */
void applyScfSafetySettings(Database::Calculation& calculation, const Utils::Settings& settings);

} // namespace DatabaseJobSubmissionHelper
} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_DATABASEJOBSUBMISSIONHELPER_H
