/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef TITRATION_DATAMANAGER_H
#define TITRATION_DATAMANAGER_H

#include <Utils/Geometry/AtomCollection.h>
#include <list>
#include <vector>

namespace Scine {
namespace StructurePreparation {

/**
 * @struct StructurePreparationData StructurePreparationData.h
 * @brief This struct holds all objects used inside the MM parametrization algorithm.
 */
struct TitrableSite {
  std::string residueName;
  int index;
  Utils::AtomCollection atoms;
  std::vector<int> indicesInFullStructure;
  bool isAcid = false;
  bool isBase = false;
  int criticalAtom;
  double refEnergy;
  double nonRefEnergy;
  double deltaE;
};

} // namespace StructurePreparation
} // namespace Scine

#endif // TITRATION_DATAMANAGER_H