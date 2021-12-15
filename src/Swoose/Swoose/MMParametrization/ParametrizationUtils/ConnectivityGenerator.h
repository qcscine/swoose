/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_CONNECTIVITYGENERATOR_H
#define MMPARAMETRIZATION_CONNECTIVITYGENERATOR_H

#include <list>
#include <memory>
#include <vector>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace Utils {
class Settings;
class BondOrderCollection;
} // namespace Utils

namespace MMParametrization {
struct ParametrizationData;

/**
 * @class ConnectivityGenerator ConnectivityGenerator.h
 * @brief This class handles the generation of the connectivity of the molecular system.
 *
 *        Right now: Bond detection with covalent radii.
 *        Possible in the future: PM6 or DFT bond orders, PDB parsing.
 */
class ConnectivityGenerator {
 public:
  /**
   * @brief Constructor.
   */
  ConnectivityGenerator(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings, Core::Log& log);
  /**
   * @brief This function generates the connectivity as a vector of lists of neighbor atom indices
   *        and updates the corresponding member in the ParametrizationData object.
   *
   *        The generated connectivity is just a guess and can be later refined after bond orders are
   *        obtained from quantum-chemical calculations of the fragments.
   */
  void generateInitialListsOfNeighbors();
  /**
   * @brief Refines the connectivity based on information obtained from quantum-chemically
   *        calculated bond orders for fragments of the whole system.
   */
  void refineListsOfNeighbors();

 private:
  // Helper function to convert a lists of neighbors vector to a bond order matrix
  Utils::BondOrderCollection generateBondOrderMatrixFromListsOfNeighbors(const std::vector<std::list<int>>& listsOfNeighbors) const;
  // The data used within all MM parametrization classes
  ParametrizationData& data_;
  // The settings
  std::shared_ptr<Utils::Settings> settings_;
  // Logger
  Core::Log& log_;
  // The threshold for a bond's bond order to be considered a bond
  double bondOrderThreshold_ = 0.4;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_CONNECTIVITYGENERATOR_H
