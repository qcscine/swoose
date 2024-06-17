/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_DATABASEORDERNAMES_H
#define MMPARAMETRIZATION_DATABASEORDERNAMES_H

namespace Scine {
namespace MMParametrization {
namespace DatabaseOrderNames {

static constexpr const char* nameOfScineBondOrdersOrder = "scine_bond_orders";
static constexpr const char* nameOfCm5ChargesOrder = "gaussian_charge_model_5";
static constexpr const char* nameOfScineHessianOrder = "scine_hessian";
static constexpr const char* nameOfScineSinglePointOrder = "scine_single_point";
static constexpr const char* nameOfScineStructureOptimizationOrder = "scine_geometry_optimization";
static constexpr const char* nameOfOrcaStructureOptimizationOrder = "orca_geometry_optimization";
static constexpr const char* nameOfTurbomoleStructureOptimizationOrder = "turbomole_geometry_optimization";

} // namespace DatabaseOrderNames
} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_DATABASEORDERNAMES_H
