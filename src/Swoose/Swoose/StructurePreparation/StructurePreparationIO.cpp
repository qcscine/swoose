/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "StructurePreparationIO.h"
#include "Protonation/TitrationData.h"
#include "SpecialCaseHandler.h"
#include "StructurePreparationData.h"
#include "StructurePreparationHelper.h"
#include <Core/Log.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/ChemicalFileFormats/FormattedStreamHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <fstream>
#include <iomanip>

namespace Scine {
namespace StructurePreparation {
namespace StructurePreparationIO {

void xyzToPdb(const std::string& xyzFile, const std::string& pdbFile) {
  Utils::AtomCollection at;
  at = Utils::ChemicalFileHandler::read(xyzFile).first;
  Utils::ChemicalFileHandler::write(pdbFile, at);
}

void writeAtomicInfoFileForProtein(const StructurePreparationData& data, const std::string& atomicInfoFile) {
  std::vector<int> listOfNegatives;
  std::vector<int> listOfPositives;

  StructurePreparationHelper::determineChargedSites(listOfNegatives, listOfPositives, data);

  if (listOfNegatives.size() == 0 && listOfPositives.size() == 0)
    return;
  std::ofstream infoFile;
  infoFile.open(atomicInfoFile);
  for (const auto& index : listOfNegatives) {
    infoFile << index << " -1 0"
             << "\n";
  }
  for (const auto& index : listOfPositives) {
    infoFile << index << " 1 0"
             << "\n";
  }
  infoFile.close();
}
void addAtomicInformationForNonRegContainer(StructurePreparationFiles& files, std::vector<std::vector<int>> subsystemMapping) {
  std::ofstream infoFile;
  infoFile.open(files.atomicInfoFile, std::ios_base::app);

  if (boost::filesystem::exists(files.nonRegContainerInfoFile) && !boost::filesystem::is_empty(files.nonRegContainerInfoFile)) {
    std::ifstream file(files.nonRegContainerInfoFile);
    for (std::string line; std::getline(file, line);) {
      std::istringstream iss(line);
      std::vector<std::string> result;
      for (std::string line; iss >> line;) {
        result.push_back(line);
      }
      if (result.size() != 3)
        throw Utils::FormattedStreamHandler::FormatMismatchException();
      int atomIndex = std::stoi(result.at(0));
      infoFile << subsystemMapping[1][atomIndex] << " " << result.at(1) << " " << result.at(2) << "\n";
    }
  }
  infoFile.close();
}

std::string getSuffix(const bfs::path& filepath) {
  std::string suffix = filepath.extension().string();

  if (suffix.size() <= 1) {
    throw Utils::FormattedStreamHandler::FormatUnsupportedException();
  }
  // Remove the leading dot
  return suffix.substr(1);
}

// TODO: This is already supported through the utils and could be removed in the future.
void writePdbFileWithResidueSpecifier(const StructurePreparationData& data, const std::string& proteinFile, Core::Log& log) {
  std::ofstream pdbFile;
  if (!data.protein.empty()) {
    pdbFile.open(proteinFile);
    assert(data.vectorOfProteinIndices.size() == data.protein.size());
    for (long unsigned int i = 0; i < data.vectorOfProteinIndices.size(); ++i) {
      int atomNumber = data.protein[i].index + 1;
      std::string residueName = data.protein[i].residueName;
      std::string atomType = data.protein[i].atomType;
      Utils::Position position = data.protein[i].position;
      Utils::ElementType element = data.fullStructure.getElement(atomNumber - 1);

      // Spaces are assigned according to the official pdb layout
      pdbFile << "ATOM" << std::setw(7) << std::right << atomNumber << "  " << std::setw(4) << std::left << atomType
              << std::setw(4) << std::left << residueName << std::setw(17) << std::right << std::fixed
              << std::setprecision(3) << position(0) << std::setw(8) << std::right << std::fixed << std::setprecision(3)
              << position(1) << std::setw(8) << std::right << std::fixed << std::setprecision(3) << position(2)
              << std::setw(7) << std::right << "    "  // occupancy
              << std::setw(6) << std::right << "     " // B factor
              << std::setw(11) << std::right << Utils::ElementInfo::symbol(element) << "\n";
    }
  }
  else
    log.output << "No protein-type substructures in your system. This structure is treated as nonRegContainer-only. "
               << Core::Log::endl;
}

} // namespace StructurePreparationIO
} // namespace StructurePreparation
} // namespace Scine
