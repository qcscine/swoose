/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AtomicInformationReader.h"
#include <Core/Log.h>
#include <regex>

namespace Scine {
namespace SwooseUtilities {

AtomicInformationReader::AtomicInformationReader(Core::Log& log) : log_(log) {
}

void AtomicInformationReader::read(const std::string& filename, std::map<int, int>& formalCharges,
                                   std::map<int, int>& unpairedElectrons, int numberOfAtoms) {
  std::ifstream indata(filename);
  if (!indata.is_open()) {
    throw std::runtime_error("The atomic information file " + filename + " cannot be opened.");
  }

  std::string line;
  while (line.empty()) {
    std::getline(indata, line);
  }

  while (!line.empty()) {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
    std::sregex_token_iterator end;

    if (*iter == "")
      iter++;

    int index = std::stoi(*iter++);
    int formalCharge = std::stoi(*iter++);
    int unpairedEle = std::stoi(*iter++);

    // Check that the number of unpaired electrons is not negative
    if (unpairedEle < 0)
      throw std::runtime_error("The number of unpaired electrons of atom " + std::to_string(index) + " cannot be negative.");

    // Check that an index was not specified more than once
    if ((formalCharges.find(index) != formalCharges.end()) || (unpairedElectrons.find(index) != unpairedElectrons.end()))
      throw std::runtime_error("In the atomic information file, don't give information about the same atom twice.");

    // Check that the index is valid
    if ((index < 0) || (index >= numberOfAtoms))
      throw std::runtime_error("At least one of the specified indices in the atomic information file is not valid.");

    if ((std::abs(formalCharge) > 5) || (unpairedEle > 6))
      log_.warning << "The formal atomic charge or the number of unpaired electrons specified for atom " << index
                   << " is very large." << Core::Log::endl;

    // Update ParametrizationData object
    if (formalCharge != 0)
      formalCharges[index] = formalCharge;
    if (unpairedEle != 0)
      unpairedElectrons[index] = unpairedEle;

    std::getline(indata, line);
  }
}

std::pair<int, int> AtomicInformationReader::read(const std::string& filename, int numberOfAtoms) {
  int totalCharge = 0;
  int totalUnpairedElectrons = 0;
  std::ifstream indata(filename);
  if (!indata.is_open()) {
    throw std::runtime_error("The atomic information file " + filename + " cannot be opened.");
  }

  std::string line;
  while (line.empty()) {
    std::getline(indata, line);
  }

  while (!line.empty()) {
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
    std::sregex_token_iterator end;

    if (*iter == "")
      iter++;

    int index = std::stoi(*iter++);
    int formalCharge = std::stoi(*iter++);
    int unpairedEle = std::stoi(*iter++);

    // Check that the number of unpaired electrons is not negative
    if (unpairedEle < 0)
      throw std::runtime_error("The number of unpaired electrons of atom " + std::to_string(index) + " cannot be negative.");

    // Check that the index is valid
    if ((index < 0) || (index >= numberOfAtoms))
      throw std::runtime_error("At least one of the specified indices in the atomic information file is not valid.");

    if ((std::abs(formalCharge) > 5) || (unpairedEle > 6))
      log_.warning << "The formal atomic charge or the number of unpaired electrons specified for atom " << index
                   << " is very large." << Core::Log::endl;

    if (formalCharge != 0)
      totalCharge += formalCharge;
    if (unpairedEle != 0)
      totalUnpairedElectrons += unpairedEle;

    std::getline(indata, line);
  }
  return std::make_pair(totalCharge, totalUnpairedElectrons);
}

} // namespace SwooseUtilities
} // namespace Scine
