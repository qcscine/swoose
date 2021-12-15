/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_ATOMICINFORMATIONREADER_H
#define MMPARAMETRIZATION_ATOMICINFORMATIONREADER_H

#include <map>
#include <string>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace SwooseUtilities {

class AtomicInformationReader {
 public:
  /**
   * @brief Constructor.
   */
  explicit AtomicInformationReader(Core::Log& log);
  /**
   * @brief Reads the formal charges and number of unpaired electrons information from file and
   *        updates the ParametrizationData object accordingly.
   * @param filename Path to the file with the information.
   * @param formalCharges Information about formal charges to update in this function.
   * @param unpairedElectrons Information about unpaired Electrons to update in this function.
   * @param numberOfAtoms Number of atoms in the system.
   */
  void read(const std::string& filename, std::map<int, int>& formalCharges, std::map<int, int>& unpairedElectrons,
            int numberOfAtoms);

 private:
  // The logger.
  Core::Log& log_;
};

} // namespace SwooseUtilities
} // namespace Scine

#endif // MMPARAMETRIZATION_ATOMICINFORMATIONREADER_H
