/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_TITRATIONFILEHANDLER_H
#define SWOOSEUTILITIES_TITRATIONFILEHANDLER_H

#include <map>
#include <string>
#include <vector>

namespace Scine {
namespace StructurePreparation {
struct TitrableSite;
} // namespace StructurePreparation

namespace MMParametrization {
struct TitrationResults;
}
namespace SwooseUtilities {
/**
 * @class ConnectivityFileHandler ConnectivityFileHandler.h
 * @brief Method for reading and writing the connectivity of atoms (list of neighbors for each atom).
 */
class TitrationFileHandler {
 public:
  TitrationFileHandler() = delete;
  /**
   * @brief
   */
  static void writeTitrationSitesFile(std::string titrationSiteFile,
                                      std::vector<StructurePreparation::TitrableSite> titrableSites);
  /**
   * @brief
   */
  static void readTitrationSitesFromFile(MMParametrization::TitrationResults& results, std::string titrationSiteFile,
                                         int numberOfAtoms, int numberOfFragments, std::map<int, std::string>& sites,
                                         std::vector<bool>& siteIsPhSensitive);
};

} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_TITRATIONFILEHANDLER_H