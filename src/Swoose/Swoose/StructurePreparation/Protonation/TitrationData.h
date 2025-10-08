/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef TITRATION_DATAMANAGER_H
#define TITRATION_DATAMANAGER_H

#include <Utils/Geometry/AtomCollection.h>
#include <vector>

namespace Scine {
namespace StructurePreparation {

/**
 * @struct StructurePreparationData StructurePreparationData.h
 * @brief This struct holds all objects used inside the MM parametrization algorithm.
 */
struct TitrableSite {
  std::string residueName;
  int index = std::numeric_limits<int>::infinity();
  Utils::AtomCollection atoms;
  std::vector<int> indicesInFullStructure;
  bool isAcid = false;
  bool isBase = false;
  int criticalAtom = std::numeric_limits<int>::infinity();
  double refEnergy = std::numeric_limits<double>::infinity();
  double nonRefEnergy = std::numeric_limits<double>::infinity();
  double deltaE = std::numeric_limits<double>::infinity();
};

} // namespace StructurePreparation
} // namespace Scine

#endif // TITRATION_DATAMANAGER_H