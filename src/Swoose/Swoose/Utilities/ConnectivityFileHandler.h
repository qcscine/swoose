/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_CONNECTIVITYFILEHANDLER_H
#define SWOOSEUTILITIES_CONNECTIVITYFILEHANDLER_H

#include <list>
#include <string>
#include <vector>

namespace Scine {
namespace SwooseUtilities {
/**
 * @class ConnectivityFileHandler ConnectivityFileHandler.h
 * @brief Method for reading and writing the connectivity of atoms (list of neighbors for each atom).
 */
class ConnectivityFileHandler {
 public:
  ConnectivityFileHandler() = delete;
  /**
   * @brief Reads the lists of neighbors from a file.
   * @param filename Name of the file from which the connectivity is read.
   * @return The connectivity of the system. It is a vector of a list of neighbor atom indices for each atom.
   */
  static std::vector<std::list<int>> readListsOfNeighbors(const std::string& filename);
  /**
   * @brief Writes all the neighbors of every atom to a file, i.e. the connectivity of the system.
   * @param filename Name of the file to which the connectivity is written.
   * @param The connectivity of the system. It is a vector of a list of neighbor atom indices for each atom.
   */
  static void writeListsOfNeighbors(const std::string& filename, const std::vector<std::list<int>>& listsOfNeighbors);
};

} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_CONNECTIVITYFILEHANDLER_H
