/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmmmCalculator.h"
#include "InteractionTermEliminator.h"
#include "QmmmCalculatorSettings.h"
#include "QmmmGradientsEvaluator.h"
#include "QmmmHelpers.h"
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Swoose/MolecularMechanics/MolecularMechanicsCalculator.h>
#include <Utils/ExternalQC/Orca/OrcaCalculatorSettings.h>
#include <iomanip>

namespace Scine {
namespace Qmmm {

std::string QmmmCalculator::name() const {
  return "QMMM";
}

bool QmmmCalculator::supportsMethodFamily(const std::string& methodFamily) const {
  return methodFamily == "QMMM";
}

QmmmCalculator::QmmmCalculator() {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<QmmmCalculatorSettings>();
}

QmmmCalculator::QmmmCalculator(const QmmmCalculator& rhs) {
  this->requiredProperties_ = rhs.requiredProperties_;
  this->setLog(rhs.getLog());
  if (this->qmCalculator_ && this->mmCalculator_) {
    this->qmCalculator_ = rhs.qmCalculator_->clone();
    this->mmCalculator_ = rhs.mmCalculator_->clone();
  }

  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));

  this->results() = rhs.results();
  this->setStructure(rhs.structure_);
  applySettings();
}

void QmmmCalculator::setUnderlyingCalculators(std::shared_ptr<Core::Calculator> qmCalculator,
                                              std::shared_ptr<Core::Calculator> mmCalculator) {
  this->mmCalculator_ = std::dynamic_pointer_cast<MolecularMechanics::MolecularMechanicsCalculator>(mmCalculator);
  if (this->mmCalculator_ == nullptr)
    throw std::runtime_error("The provided MM calculator is not valid.");

  this->qmCalculator_ = qmCalculator;
  setLogForUnderlyingCalculators();
  // Pass the QM and MM calculator settings to the QM/MM settings // TODO: Should this be done?
  auto& qmmmSettings = static_cast<QmmmCalculatorSettings&>(*this->settings_);
  qmmmSettings.addExternalSettings(this->qmCalculator_->settings());
  qmmmSettings.addExternalSettings(this->mmCalculator_->settings());
}

void QmmmCalculator::applySettings() {
  using namespace SwooseUtilities::SettingsNames;
  settings_->normalizeStringCases(); // convert all option names to lower case letters
  if (settings_->valid()) {
    setLogForUnderlyingCalculators();
    ignoreQm_ = settings_->getBool(ignoreQmOption);
    qmRegionFile_ = settings_->getString(qmRegionXyzFile);
    calculateReducedQmMmEnergy_ = settings_->getBool(calculateReducedQmMmEnergy);
    chargeRedistributionScheme_ = settings_->getString(chargeRedistributionKey);
    // The electrostatic embedding setting has to be set before the calculators get their settings updated.
    electrostaticEmbedding_ = settings_->getBool(electrostaticEmbedding);
    applySettingsToUnderlyingCalculators();
    // Process the QM atoms setting:
    auto qmAtoms = settings_->getIntList(qmAtomsList);
    if (qmAtoms != listOfQmAtoms_) { // Only update the vector if necessary.
      if (qmAtoms.empty())
        throw std::runtime_error("Please specify at least one atom for the QM region.");
      listOfQmAtoms_ = qmAtoms;
      QmmmHelpers::checkValidityOfQmRegion(listOfQmAtoms_, structure_);
    }
  }
  else {
    settings_->throwIncorrectSettings();
  }
}

const Utils::Results& QmmmCalculator::calculate(std::string description) {
  applySettings();
  try {
    return calculateImpl(description);
  }
  catch (std::runtime_error& e) {
    throw Core::UnsuccessfulCalculationException(e.what());
  }
}

void QmmmCalculator::setStructure(const Utils::AtomCollection& structure) {
  if (mmCalculator_ == nullptr || qmCalculator_ == nullptr)
    throw std::runtime_error(
        "The QM/MM calculator must be provided with a QM and an MM calculator before setting a structure.");

  structure_ = structure;
  scfConvCriterionIsSet_ = false;
  applySettings();
  mmCalculator_->setStructure(structure_);

  mmBoundaryAtoms_.clear();
  qmRegion_ = QmmmHelpers::createQmRegion(listOfQmAtoms_, structure_, mmCalculator_->listsOfNeighbors(), qmRegionFile_,
                                          mmBoundaryAtoms_);
  // Only set the structure for the QM calculator if no calculation has been conducted yet and if the elements
  // of the structure changed. This prevents too many calculation directories during a geometry optimization.
  if (results_.has<Utils::Property::SuccessfulCalculation>() &&
      qmCalculator_->getStructure()->getElements() == qmRegion_.getElements()) {
    qmCalculator_->modifyPositions(qmRegion_.getPositions());
  }
  else {
    qmCalculator_->setStructure(qmRegion_);
  }

  prepareTermsForMmCalculator();

  if (electrostaticEmbedding_) {
    auto redistributedCharges =
        QmmmHelpers::getRedistributedCharges(mmCalculator_->atomicCharges(), structure_.getPositions(), mmBoundaryAtoms_,
                                             mmCalculator_->listsOfNeighbors(), listOfQmAtoms_, chargeRedistributionScheme_);
    QmmmHelpers::writePointChargesFile(structure_.getPositions(), redistributedCharges, listOfQmAtoms_, pointChargesFilename_);
  }
  else {
    boost::filesystem::remove(pointChargesFilename_);
  }
}

const Utils::Results& QmmmCalculator::calculateImpl(std::string description) {
  this->mmCalculator_->setLog(this->getLog());

  if (ignoreQm_)
    return ignoreQmCalculateImpl(description);

  if (qmCalculator_->getStructure()->size() == 0)
    throw std::runtime_error("No QM atoms were specified!");

  this->getLog().output << "\n\n";
  this->getLog().output << "Calculating QM energy..." << Core::Log::endl;
  auto qmResults = qmCalculator_->calculate("QM calculation for QM/MM");
  double qmEnergy = qmResults.get<Utils::Property::Energy>();
  this->getLog().output << "Calculating MM energy..." << Core::Log::endl;
  auto mmResults = mmCalculator_->calculate("MM calculation for QM/MM");
  double mmEnergy = mmResults.get<Utils::Property::Energy>();
  auto mmGradients = mmResults.get<Utils::Property::Gradients>();

  // Total energy:
  double totalEnergy = qmEnergy + mmEnergy;

  // Output:
  this->getLog().debug << std::setprecision(10) << "QM energy (Hartree): " << qmEnergy << Core::Log::endl;
  this->getLog().debug << "MM energy (Hartree): " << mmEnergy << Core::Log::endl;

  if (calculateReducedQmMmEnergy_) {
    prepareTermsForMmCalculator(true);
    mmResults = mmCalculator_->calculate("MM calculation for QM/MM, no environment-only contributions");
    double reducedEnergy = qmEnergy + mmResults.get<Utils::Property::Energy>();
    this->getLog().output << std::setprecision(10) << "Reduced QM/MM energy (Hartree): " << reducedEnergy << Core::Log::endl;
    prepareTermsForMmCalculator(false);
  }

  // Assemble results
  results_.set<Utils::Property::Description>(std::move(description));
  if (requiredProperties_.containsSubSet(Utils::Property::Energy))
    results_.set<Utils::Property::Energy>(totalEnergy);
  if (requiredProperties_.containsSubSet(Utils::Property::Gradients)) {
    Utils::GradientCollection pcGradientsContributions;
    if (electrostaticEmbedding_)
      pcGradientsContributions = qmResults.get<Utils::Property::PointChargesGradients>();

    QmmmGradientsEvaluator gradEvaluator(qmResults.get<Utils::Property::Gradients>(), mmGradients,
                                         pcGradientsContributions, listOfQmAtoms_, mmBoundaryAtoms_,
                                         mmCalculator_->listsOfNeighbors(), structure_, qmRegion_);
    results_.set<Utils::Property::Gradients>(gradEvaluator.calculateQmmmGradients());
  }

  results_.set<Utils::Property::SuccessfulCalculation>(true);

  return results_;
}

std::unique_ptr<Utils::AtomCollection> QmmmCalculator::getStructure() const {
  return std::make_unique<Utils::AtomCollection>(structure_);
}

void QmmmCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  structure_.setPositions(newPositions);
  mmCalculator_->modifyPositions(std::move(newPositions));

  mmBoundaryAtoms_.clear();
  qmRegion_ = QmmmHelpers::createQmRegion(listOfQmAtoms_, structure_, mmCalculator_->listsOfNeighbors(), qmRegionFile_,
                                          mmBoundaryAtoms_);
  qmCalculator_->modifyPositions(qmRegion_.getPositions());

  if (electrostaticEmbedding_) {
    auto redistributedCharges =
        QmmmHelpers::getRedistributedCharges(mmCalculator_->atomicCharges(), structure_.getPositions(), mmBoundaryAtoms_,
                                             mmCalculator_->listsOfNeighbors(), listOfQmAtoms_, chargeRedistributionScheme_);
    QmmmHelpers::writePointChargesFile(structure_.getPositions(), redistributedCharges, listOfQmAtoms_, pointChargesFilename_);
  }
  else {
    boost::filesystem::remove(pointChargesFilename_);
  }
}

const Utils::PositionCollection& QmmmCalculator::getPositions() const {
  return structure_.getPositions();
}

void QmmmCalculator::setRequiredProperties(const Utils::PropertyList& requiredProperties) {
  requiredProperties_ = requiredProperties;
  mmCalculator_->setRequiredProperties(requiredProperties);

  // This setting has to be recopied from the settings object, because it might have changed recently:
  electrostaticEmbedding_ = settings_->getBool(SwooseUtilities::SettingsNames::electrostaticEmbedding);

  if (qmCalculator_) {
    if (electrostaticEmbedding_ && requiredProperties.containsSubSet(Utils::Property::Gradients)) {
      auto rP = requiredProperties;
      rP.addProperty(Utils::Property::PointChargesGradients);
      qmCalculator_->setRequiredProperties(rP);
    }
    else {
      qmCalculator_->setRequiredProperties(requiredProperties);
    }
  }
}

Utils::PropertyList QmmmCalculator::getRequiredProperties() const {
  return requiredProperties_;
}

Utils::PropertyList QmmmCalculator::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients;
}

const Utils::Settings& QmmmCalculator::settings() const {
  return *settings_;
}

Utils::Settings& QmmmCalculator::settings() {
  return *settings_;
}

const Utils::Results& QmmmCalculator::results() const {
  return results_;
}

Utils::Results& QmmmCalculator::results() {
  return results_;
}

void QmmmCalculator::loadState(std::shared_ptr<Core::State> /*state*/) {
  // empty implementation
}

std::shared_ptr<Core::State> QmmmCalculator::getState() const {
  return nullptr;
}

void QmmmCalculator::applySettingsToUnderlyingCalculators() {
  // MM calculator
  if (mmCalculator_) {
    for (const auto& descriptor : mmCalculator_->settings().getDescriptorCollection()) {
      auto key = descriptor.first;
      if (settings_->valueExists(key))
        mmCalculator_->settings().modifyValue(key, settings_->getValue(key));
    }
  }

  // QM calculator
  if (qmCalculator_) {
    if (electrostaticEmbedding_) {
      if (qmCalculator_->settings().valueExists(Utils::ExternalQC::SettingsNames::pointChargesFile))
        qmCalculator_->settings().modifyString(Utils::ExternalQC::SettingsNames::pointChargesFile, pointChargesFilename_);
      else
        throw std::runtime_error("The selected calculator does not support electrostatic embedding.");
    }
    for (const auto& d : qmCalculator_->settings().getDescriptorCollection()) {
      if (!this->settings().valueExists(d.first))
        continue;
      if (d.first == Utils::ExternalQC::SettingsNames::pointChargesFile)
        continue;
      if (d.first == Utils::SettingsNames::selfConsistenceCriterion) {
        if (scfConvCriterionIsSet_)
          continue;
        scfConvCriterionIsSet_ = true;
      }
      qmCalculator_->settings().modifyValue(d.first, settings_->getValue(d.first));
    }
  }
}

void QmmmCalculator::prepareTermsForMmCalculator(bool reducedEnergyCalculation) {
  InteractionTermEliminator termsEliminator(listOfQmAtoms_, mmCalculator_);
  termsEliminator.reset();
  termsEliminator.eliminateInteractionTerms(electrostaticEmbedding_, reducedEnergyCalculation);
}

void QmmmCalculator::setLogForUnderlyingCalculators() {
  if (this->mmCalculator_ != nullptr)
    this->mmCalculator_->setLog(this->getLog());
  if (this->qmCalculator_ != nullptr) {
    Core::Log warningLog = Core::Log::silent();
    warningLog.warning.add("cerr", Core::Log::cerrSink());
    warningLog.error.add("cerr", Core::Log::cerrSink());
    this->qmCalculator_->setLog(warningLog);
  }
}

QmmmCalculator::~QmmmCalculator() {
  boost::filesystem::remove(pointChargesFilename_);
}

const Utils::Results& QmmmCalculator::ignoreQmCalculateImpl(std::string description) {
  this->getLog().output << "\n\n";
  this->getLog().output << "Ignoring QM calculation..." << Core::Log::nl << "Calculating MM energy..." << Core::Log::endl;
  auto mmResults = mmCalculator_->calculate("MM calculation for QM/MM");

  // Assemble results
  results_.set<Utils::Property::Description>(std::move(description));
  if (requiredProperties_.containsSubSet(Utils::Property::Energy))
    results_.set<Utils::Property::Energy>(mmResults.get<Utils::Property::Energy>());
  if (requiredProperties_.containsSubSet(Utils::Property::Gradients)) {
    results_.set<Utils::Property::Gradients>(mmResults.get<Utils::Property::Gradients>());
  }

  results_.set<Utils::Property::SuccessfulCalculation>(true);

  return results_;
}

} // namespace Qmmm
} // namespace Scine
