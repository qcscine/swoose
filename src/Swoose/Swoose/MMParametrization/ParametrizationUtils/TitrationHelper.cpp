/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TitrationHelper.h"
#include "Swoose/StructurePreparation/Protonation/ProtonationHandler.h"
#include "Swoose/StructurePreparation/Protonation/ProtonationHelper.h"
#include <Swoose/MMParametrization/ParametrizationData.h>
#include <Swoose/Utilities/AtomicInformationReader.h>
#include <Swoose/Utilities/SettingsNames.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/NativeFilenames.h>
#include <boost/filesystem.hpp>

namespace Scine {
namespace MMParametrization {
using namespace StructurePreparation;

TitrationHelper::TitrationHelper(std::shared_ptr<Utils::Settings>& settings) : settings_(settings) {
  // todo
}

Utils::AtomCollection TitrationHelper::changeProtonationState(const Utils::AtomCollection& refStructure,
                                                              std::string residueName, bool isBase, int indexOfCriticalAtom,
                                                              std::vector<int> superfluousHydrogens) {
  Utils::AtomCollection newStructure;
  newStructure = refStructure;
  ProtonationHelper::removeProtonsFromStructure(newStructure, superfluousHydrogens);
  if (isBase) {
    ProtonationHandler handler;
    handler.reprotonateCriticalAtom(newStructure, indexOfCriticalAtom, residueName);
  }
  return newStructure;
}

double TitrationHelper::getModelPka(const std::string& residueName) {
  AminoAcidCategorizer categorizer;
  auto it = categorizer.modelPkaMap.find(residueName);
  if (it != categorizer.modelPkaMap.end())
    return it->second;
  throw std::runtime_error("model pka to " + residueName + " could not be determined");
}

void TitrationHelper::collectTrainingData(TitrationResults& results) {
  for (const auto& availableGroup : results.requiredFunctionalGroups) {
    TrainingData data;
    data.groupType = availableGroup;
    auto trainingDataDirectory = settings_->getString(SwooseUtilities::SettingsNames::trainingDataDirectory);
    if (!boost::filesystem::exists(Utils::NativeFilenames::combinePathSegments(trainingDataDirectory, availableGroup)))
      throw std::runtime_error("Invalid training data directory! Cannot find training data for " + availableGroup + " group.");
    boost::filesystem::directory_iterator iterator(
        Utils::NativeFilenames::combinePathSegments(trainingDataDirectory, availableGroup));
    while (iterator != boost::filesystem::directory_iterator{}) {
      std::string directory = iterator->path().string();
      data.pKas.push_back(getPkaOfTrainingMolecule(directory));
      data.eneryOfDeprotonation.push_back(getEnergyOfDeprotonationForTrainingMolecule(directory));
      ++iterator;
    }
    results.trainingData.push_back(data);
  }
}

void TitrationHelper::calculateFreeEnergiesOfDeprotonation(TitrationResults& results) {
  for (auto& site : results.sites) {
    if (site.isBase)
      site.deltaE = site.refEnergy - site.nonRefEnergy;
    else if (site.isAcid)
      site.deltaE = site.nonRefEnergy - site.refEnergy;
  }
}

double TitrationHelper::getPkaOfTrainingMolecule(std::string dataDirectory) {
  std::string pkaFile = Utils::NativeFilenames::combinePathSegments(dataDirectory, "pka.dat");
  std::ifstream indata(pkaFile);
  if (!indata.is_open()) {
    throw std::runtime_error("The pka file " + pkaFile + " cannot be opened.");
  }

  std::string line;
  double pka = 0.0;
  try {
    std::getline(indata, line);
    pka = std::stod(line);
  }
  catch (std::exception& e) {
    throw std::runtime_error("No pka given for training fragment in " + dataDirectory);
  }
  return pka;
}

double TitrationHelper::getEnergyOfDeprotonationForTrainingMolecule(std::string dataDirectory) {
  (void)dataDirectory;
  return 0.0;
}

} // namespace MMParametrization
} // namespace Scine
