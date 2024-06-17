/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TitrationFileHandler.h"
#include <Swoose/MMParametrization/ParametrizationData.h>
#include <Swoose/StructurePreparation/ProteinStructures.h>
#include <Swoose/StructurePreparation/Protonation/TitrationData.h>
#include <fstream>
#include <regex>
#include <string>

namespace Scine {
namespace SwooseUtilities {

void TitrationFileHandler::readTitrationSitesFromFile(MMParametrization::TitrationResults& results,
                                                      std::string titrationSiteFile, int numberOfAtoms,
                                                      int numberOfFragments, std::map<int, std::string>& sites,
                                                      std::vector<bool>& siteIsPhSensitive) {
  std::ifstream indata(titrationSiteFile);
  if (!indata.is_open()) {
    throw std::runtime_error("The titration site file " + titrationSiteFile + " cannot be opened.");
  }

  std::string line;
  std::regex r1("([0-9]+) ([a-zA-Z]+)");
  std::smatch m1;
  while (std::getline(indata, line)) {
    // check if format is correct
    if (std::regex_search(line, m1, r1)) {
      StructurePreparation::TitrableSite site;
      site.residueName = m1[2];
      site.index = std::stoi(m1[1]);
      // Check that the index is valid
      if ((site.index < 0) || (site.index >= numberOfAtoms))
        throw std::runtime_error("At least one of the specified indices in the titration sites file is not valid.");
      StructurePreparation::AminoAcidCategorizer categorizer;

      if (std::find(categorizer.acids.begin(), categorizer.acids.end(), site.residueName) != categorizer.acids.end())
        site.isAcid = true;
      else if (std::find(categorizer.bases.begin(), categorizer.bases.end(), site.residueName) != categorizer.bases.end())
        site.isBase = true;
      if (!site.isAcid && !site.isBase)
        throw std::runtime_error("The site " + site.residueName + " is none of the known titrable amino acids. ");

      auto functionalGroup = categorizer.functionalGroups.find(site.residueName)->second;
      // save the type of functional group for regression
      results.requiredFunctionalGroups.push_back(functionalGroup);
      // save the site
      results.sites.push_back(site);
      sites.insert({site.index, site.residueName});
      if (numberOfFragments == 1)
        siteIsPhSensitive.at(0) = true;
      else
        siteIsPhSensitive.at(site.index) = true;
    }
    else {
      throw std::runtime_error("Incorrect file format!");
    }
  }
}

void TitrationFileHandler::writeTitrationSitesFile(std::string titrationSiteFile,
                                                   std::vector<StructurePreparation::TitrableSite> titrableSites) {
  std::ofstream file;
  file.open(titrationSiteFile);

  for (auto& site : titrableSites) {
    file << site.criticalAtom << " " << site.residueName << "\n";
  }
  file.close();
}

} // namespace SwooseUtilities
} // namespace Scine