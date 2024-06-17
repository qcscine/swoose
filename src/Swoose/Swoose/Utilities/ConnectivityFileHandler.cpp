/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ConnectivityFileHandler.h"
#include <fstream>
#include <regex>

namespace Scine {
namespace SwooseUtilities {

std::vector<std::list<int>> ConnectivityFileHandler::readListsOfNeighbors(const std::string& filename) {
  std::vector<std::list<int>> listsOfNeighbors;

  std::ifstream indata(filename);
  if (!indata.is_open()) {
    throw std::runtime_error("The connectivity file " + filename + " cannot be opened.");
  }

  std::string line;
  while (line.empty()) {
    std::getline(indata, line);
  }

  while (!line.empty()) {
    std::list<int> listForOneAtom;
    std::regex rgx("\\s+");
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
    std::sregex_token_iterator end;

    if (*iter == "")
      iter++;

    while (iter != end) {
      int atomIndex = std::stoi(*iter++);
      if (int(listsOfNeighbors.size()) == atomIndex) // atom cannot be bonded to itself
        throw std::runtime_error("Error in connectivity file. Atom " + std::to_string(atomIndex) + " is bonded to itself.");
      if (atomIndex >= 0)
        listForOneAtom.push_back(atomIndex);
    }

    listsOfNeighbors.push_back(listForOneAtom);
    if (indata.eof())
      break;
    std::getline(indata, line);
  }

  // Check for validity
  for (int i = 0; i < int(listsOfNeighbors.size()); ++i) {
    for (const auto& neighbor : listsOfNeighbors[i]) {
      bool indexIsValid = neighbor >= 0 && neighbor < int(listsOfNeighbors.size());
      bool foundInCorrespondingPlace = std::find(listsOfNeighbors[neighbor].begin(), listsOfNeighbors[neighbor].end(), i) !=
                                       listsOfNeighbors[neighbor].end();
      if (!indexIsValid || !foundInCorrespondingPlace)
        throw std::runtime_error("Connectivity file is invalid! Error during check of atom with index: " + std::to_string(i));
    }
  }

  return listsOfNeighbors;
}

void ConnectivityFileHandler::writeListsOfNeighbors(const std::string& filename,
                                                    const std::vector<std::list<int>>& listsOfNeighbors) {
  std::ofstream connectivityFile(filename);

  for (int i = 0; i < int(listsOfNeighbors.size()); ++i) {
    for (const auto& neighbor : listsOfNeighbors[i]) {
      connectivityFile << neighbor << "  ";
    }
    if (listsOfNeighbors[i].empty())
      connectivityFile << "-1  ";
    connectivityFile << std::endl;
  }
}

} // namespace SwooseUtilities
} // namespace Scine
