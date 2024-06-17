/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MolecularSystemPartitioner.h"
#include "MMParametrizationSettings.h"
#include "ParametrizationData.h"
#include <Core/Log.h>
#include <Swoose/Utilities/FragmentAnalyzer.h>
#include <Swoose/Utilities/FragmentationHelper.h>
#include <Swoose/Utilities/SubsystemGenerator.h>
#include <Utils/Constants.h>
#include <iomanip>

namespace Scine {
namespace MMParametrization {

MolecularSystemPartitioner::MolecularSystemPartitioner(ParametrizationData& data,
                                                       std::shared_ptr<Utils::Settings> settings, Core::Log& log)
  : data_(data),
    settings_(settings),
    constrainedAtomsIdentifier_(std::make_unique<ConstrainedAtomsIdentifier>(data)),
    log_(log) {
  numberAtomsThreshold_ = settings_->getInt(SwooseUtilities::SettingsNames::numberAtomsThreshold);
  subsystemRadius_ = settings_->getDouble(SwooseUtilities::SettingsNames::subsystemRadius) * Utils::Constants::bohr_per_angstrom;
  bondOrderThreshold_ = settings_->getDouble(SwooseUtilities::SettingsNames::bondOrderThreshold);
}

void MolecularSystemPartitioner::divideSystem() {
  if (data_.numberOfAtoms <= numberAtomsThreshold_)
    prepareDataForOneSubsystem();
  else
    divideSystemIntoSubsystems();
}

void MolecularSystemPartitioner::prepareDataForOneSubsystem() {
  SwooseUtilities::FragmentAnalyzer fragmentAnalyzer(data_.formalCharges, data_.unpairedElectrons);
  fragmentAnalyzer.analyzeFragment(data_.fullStructure);
  data_.vectorOfStructures.push_back(std::make_unique<Utils::AtomCollection>(data_.fullStructure));
  data_.vectorOfChargesAndMultiplicities.emplace_back(
      std::make_pair(fragmentAnalyzer.getMolecularCharge(), fragmentAnalyzer.getSpinMultiplicity()));

  // Update index mapping information
  data_.atomIndexMapping.resize(1);
  SwooseUtilities::FragmentationHelper::updateInformationForIndexMapping(data_.fullStructure, data_.fullStructure,
                                                                         data_.atomIndexMapping[0]);

  // Check validity of the whole system
  bool valid = fragmentAnalyzer.analyzeFragment(data_.fullStructure);
  if (!valid)
    throw std::runtime_error(
        "The given system with its formal charges and numbers of unpaired electrons is not valid.");
}

void MolecularSystemPartitioner::divideSystemIntoSubsystems() {
  // Construct fragment analyzer
  SwooseUtilities::FragmentAnalyzer fragmentAnalyzer(data_.formalCharges, data_.unpairedElectrons);

  // Resize the atom index mapping vector
  data_.atomIndexMapping.resize(data_.numberOfAtoms);
  // Reserve some space for all the structure pointers
  data_.vectorOfStructures.reserve(data_.numberOfAtoms);
  // Loop over all atoms in the full structure

  SwooseUtilities::SubsystemGenerator subsystemGenerator(data_.fullStructure, data_.bondOrders, fragmentAnalyzer,
                                                         bondOrderThreshold_, numberAtomsThreshold_, log_);

  for (int atomIndex = 0; atomIndex < data_.fullStructure.size(); ++atomIndex) {
    auto subsystem = subsystemGenerator.generateSubsystem(atomIndex, data_.atomIndexMapping.at(atomIndex), subsystemRadius_);
    data_.vectorOfStructures.push_back(std::make_unique<Utils::AtomCollection>(std::move(subsystem)));
    data_.vectorOfChargesAndMultiplicities.emplace_back(
        std::make_pair(fragmentAnalyzer.getMolecularCharge(), fragmentAnalyzer.getSpinMultiplicity()));
  }

  // Update information about which atoms to constrain during a geometry optimization
  int atomIndex = 0;
  for (const auto& subsystem : data_.vectorOfStructures) {
    constrainedAtomsIdentifier_->updateInformationAboutConstrainedAtoms(*subsystem, atomIndex);
    atomIndex++;
  }
}

} // namespace MMParametrization
} // namespace Scine
