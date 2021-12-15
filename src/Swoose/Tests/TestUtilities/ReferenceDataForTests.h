/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_REFERENCEDATAFORTESTS_H
#define SWOOSE_REFERENCEDATAFORTESTS_H

#include "../Files/tests_file_location.h"
#include <Eigen/Core>
#include <fstream>
#include <list>
#include <regex>
#include <vector>

namespace Scine {
namespace Swoose {
using ForcesCollection = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

namespace ReferenceDataForTests {

inline std::vector<double> parseReferenceEnergies() {
  std::vector<double> energies;

  std::ifstream indata(reference_energies_file);
  if (!indata.is_open()) {
    throw std::runtime_error("The file containing the reference energies for the tests could not be opened.");
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

    double value = std::stod(*iter);
    energies.push_back(value);

    std::getline(indata, line);
  }

  return energies;
}

inline std::vector<ForcesCollection> parseReferenceForces(std::string filename, int nStruct, int nAtoms) {
  std::vector<Swoose::ForcesCollection> forces;
  for (int i = 0; i < nStruct; ++i) {
    ForcesCollection f(nAtoms, 3);
    forces.push_back(f);
  }

  std::ifstream indata(filename);
  if (!indata.is_open()) {
    throw std::runtime_error("The file containing the reference forces for the tests could not be opened.");
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

    int structureIndex = std::stoi(*iter++);
    int atomIndex = std::stoi(*iter++);
    int dimension = std::stoi(*iter++);
    double value = std::stod(*iter++);

    forces.at(structureIndex)(atomIndex, dimension) = value;

    std::getline(indata, line);
  }

  return forces;
}

} // namespace ReferenceDataForTests
} // namespace Swoose
} // namespace Scine

#endif // SWOOSE_REFERENCEDATAFORTESTS_H
