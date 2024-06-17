/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmmmCalculator.h"
#include "InteractionTermEliminator.h"
#include "QmmmCalculatorSettings.h"
#include "QmmmGradientsEvaluator.h"
#include "QmmmHelpers.h"
#include "QmmmHessianEvaluator.h"
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Swoose/MolecularMechanics/MolecularMechanicsCalculator.h>
#include <Utils/ExternalQC/Orca/OrcaCalculatorSettings.h>
#include <Utils/GeometryOptimization/CoordinateSystem.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/Optimizer/GradientBased/Bfgs.h>
#include <Utils/UniversalSettings/OptimizationSettingsNames.h>
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

QmmmCalculator::QmmmCalculator(const QmmmCalculator& rhs)
  : CloneInterface(rhs), requiredProperties_(rhs.requiredProperties_) {
  this->setLog(rhs.getLog());
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  setUnderlyingCalculatorsImpl(rhs.getUnderlyingCalculators());
  valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  setStructureImpl(rhs.structure_);
  applySettings();
  results_ = rhs.results_; // after structure is set
}

void QmmmCalculator::setUnderlyingCalculators(std::vector<std::shared_ptr<Core::Calculator>> underlyingCalculators) {
  setUnderlyingCalculatorsImpl(underlyingCalculators);
}

void QmmmCalculator::setUnderlyingCalculatorsImpl(std::vector<std::shared_ptr<Core::Calculator>> underlyingCalculators) {
  if (underlyingCalculators.size() != 2) {
    throw std::runtime_error("The QM/MM calculator needs exactly two underlying calculators; received " +
                             std::to_string(underlyingCalculators.size()) + " instead.");
  }
  this->mmCalculator_ =
      std::dynamic_pointer_cast<MolecularMechanics::MolecularMechanicsCalculator>(underlyingCalculators.at(1)->clone());
  if (this->mmCalculator_ == nullptr)
    throw std::runtime_error("The provided MM calculator is not valid.");

  this->qmCalculator_ = underlyingCalculators.at(0)->clone();
  if (this->qmCalculator_ == nullptr)
    throw std::runtime_error("The provided QM calculator is not valid.");
  // get rid of all Calculator-specific settings
  removeCalculatorSpecificSettings();
  // add the settings of the new calculators
  addUnderlyingSettingsImpl();
}

void QmmmCalculator::addUnderlyingSettings() {
  addUnderlyingSettingsImpl();
}

void QmmmCalculator::addUnderlyingSettingsImpl() {
  // Pass the QM and MM calculator settings to the QM/MM settings
  auto& qmmmSettings = dynamic_cast<QmmmCalculatorSettings&>(*this->settings_);
  if (this->mmCalculator_) {
    qmmmSettings.addExternalSettings(this->mmCalculator_->settings());
  }
  if (this->qmCalculator_) {
    qmmmSettings.addExternalSettings(this->qmCalculator_->settings());
  }
}

void QmmmCalculator::removeCalculatorSpecificSettings() {
  const auto previousSettings = Utils::Settings(*settings_);
  this->settings_ = std::make_unique<QmmmCalculatorSettings>();
  this->settings_->merge(previousSettings, true); // merge all duplicate keys (only QM/MM base settings)
}

std::vector<std::shared_ptr<Core::Calculator>> QmmmCalculator::getUnderlyingCalculators() const {
  return std::vector<std::shared_ptr<Core::Calculator>>{qmCalculator_, mmCalculator_};
}

void QmmmCalculator::applySettings() {
  using namespace SwooseUtilities::SettingsNames;
  settings_->normalizeStringCases(); // convert all option names to lower case letters
  if (settings_->valid()) {
    ignoreQm_ = settings_->getBool(Utils::SettingsNames::ignoreQmOption);
    qmRegionFile_ = settings_->getString(qmRegionXyzFile);
    calculateReducedQmMmEnergy_ = settings_->getBool(calculateReducedQmMmEnergy);
    chargeRedistributionScheme_ = settings_->getString(chargeRedistributionKey);
    // The electrostatic embedding setting has to be set before the calculators get their settings updated.
    electrostaticEmbedding_ = settings_->getBool(Utils::SettingsNames::electrostaticEmbedding);
    turbomoleIsQmCalculator_ = (qmCalculator_->name() == "TURBOMOLE");
    optimizeLinks_ = settings_->getBool(Utils::SettingsNames::optimizeLinks);
    if (ignoreQm_ && optimizeLinks_) {
      getLog().warning << "The QM region is ignored, the option to optimize QM link atoms has no effect" << Core::Log::nl;
    }
    applySettingsToUnderlyingCalculators();
    // Process the QM atoms setting:
    auto qmAtoms = settings_->getIntList(Utils::SettingsNames::qmAtomsList);
    if (qmAtoms != listOfQmAtoms_) { // Only update the vector if necessary.
      if (qmAtoms.empty())
        throw std::runtime_error("Please specify at least one atom for the QM region.");
      listOfQmAtoms_ = qmAtoms;
      QmmmHelpers::checkValidityOfQmRegion(listOfQmAtoms_, structure_);
    }
    mmAtomsLeft_ = size_t(this->structure_.size()) > this->listOfQmAtoms_.size();
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
  setStructureImpl(structure);
}

void QmmmCalculator::setStructureImpl(const Utils::AtomCollection& structure) {
  if (mmCalculator_ == nullptr)
    throw std::runtime_error("The QM/MM calculator must be provided with an MM calculator before setting a structure.");
  if (qmCalculator_ == nullptr)
    throw std::runtime_error(
        "The QM/MM calculator must be provided with a valid QM calculator before setting a structure.");
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

  handleElectrostaticEmbedding();
}

void QmmmCalculator::handleElectrostaticEmbedding() {
  if (!mmCalculator_) {
    return;
  }
  if (electrostaticEmbedding_) {
    auto redistributedCharges =
        QmmmHelpers::getRedistributedCharges(mmCalculator_->atomicCharges(), structure_.getPositions(), mmBoundaryAtoms_,
                                             mmCalculator_->listsOfNeighbors(), listOfQmAtoms_, chargeRedistributionScheme_);
    if (qmCalculator_->settings().valueExists(Utils::SettingsNames::mmCharges)) {
      if (settings().getString(SwooseUtilities::SettingsNames::chargeRedistributionKey) ==
          SwooseUtilities::OptionNames::redistributedChargeAndDipolesOption) {
        std::string error = "The '" + std::string(Utils::SettingsNames::mmCharges) +
                            "' option is not compatible with the charge redistribution scheme '" +
                            std::string{SwooseUtilities::OptionNames::redistributedChargeAndDipolesOption} + "'.";
        if (!qmCalculator_->settings().valueExists(Utils::ExternalQC::SettingsNames::pointChargesFile)) {
          // we have no other chance to set electrostatics, so we error out
          throw std::runtime_error(error);
        }
        // we should have another chance, so we simply warn about it
        this->getLog().warning << error << Core::Log::endl;
      }
      else {
        // we should be able to directly set the charges
        const auto chargesAndPositions =
            QmmmHelpers::writeChargesAndPositionsAsList(structure_, redistributedCharges, listOfQmAtoms_);
        qmCalculator_->settings().modifyDoubleList(Utils::SettingsNames::mmCharges, chargesAndPositions);
      }
    }
    else if (qmCalculator_->settings().valueExists(Utils::ExternalQC::SettingsNames::pointChargesFile)) {
      // the turbomole calculator requires a different point charges file format
      QmmmHelpers::writePointChargesFile(structure_.getPositions(), redistributedCharges, listOfQmAtoms_,
                                         pointChargesFilename_, turbomoleIsQmCalculator_);
      qmCalculator_->settings().modifyString(Utils::ExternalQC::SettingsNames::pointChargesFile, pointChargesFilename_);
    }
    else {
      this->getLog().warning
          << "The selected calculator does not support electrostatic embedding. Switching off electrostatic embedding"
          << Core::Log::endl;
      settings_->modifyBool(Utils::SettingsNames::electrostaticEmbedding, false);
      electrostaticEmbedding_ = false;
    }
  }
  if (!electrostaticEmbedding_) {
    boost::filesystem::remove(pointChargesFilename_);
  }
}

void QmmmCalculator::optimizeLinks() {
  if (listOfQmAtoms_.size() == size_t(qmCalculator_->getStructure()->size())) {
    getLog().debug << "There are no link atoms, no optimization necessary" << Core::Log::endl;
    return;
  }
  auto requiredProperties = qmCalculator_->getRequiredProperties();
  Utils::CalculationRoutines::setLog(*qmCalculator_, true, true, true);
  auto optimizer = Utils::GeometryOptimizer<Utils::Bfgs>(*qmCalculator_);
  int maxCycles = 50;
  optimizer.check.maxIter = maxCycles;
  optimizer.check.deltaValue = 1e-6;
  std::vector<int> fixedAtoms;
  for (int i = 0; i < static_cast<int>(listOfQmAtoms_.size()); ++i) {
    fixedAtoms.push_back(i);
  }
  // add observer to have trajectory
  std::ofstream trajectory("link_optimization.trj.xyz", std::ofstream::out);
  auto atoms = *(qmCalculator_->getStructure());
  auto func = [&](const int& /* cycle */, const double& energy, const Eigen::VectorXd& /* params */) {
    Utils::XyzStreamHandler::write(trajectory, atoms, std::to_string(energy));
  };
  optimizer.addObserver(func);

  // this works because link atoms are always at the end of the list
  optimizer.fixedAtoms = fixedAtoms;
  optimizer.coordinateSystem = Utils::CoordinateSystem::Cartesian;
  int nCycles = 0;
  try {
    nCycles = optimizer.optimize(atoms, getLog());
  }
  catch (...) {
    throw;
  }
  if (nCycles == maxCycles) {
    getLog().warning << "QM/MM link optimization did not converge within " << maxCycles << " cycles." << Core::Log::endl;
  }
  qmCalculator_->modifyPositions(atoms.getPositions());
  qmCalculator_->setRequiredProperties(requiredProperties);
}

const Utils::Results& QmmmCalculator::calculateImpl(std::string description) {
  setLogForUnderlyingCalculators();

  if (ignoreQm_)
    return ignoreQmCalculateImpl(description);

  if (qmCalculator_->getStructure()->size() == 0)
    throw std::runtime_error("No QM atoms were specified!");

  if (optimizeLinks_) {
    this->getLog().debug << "Optimizing QM links..." << Core::Log::endl;
    optimizeLinks();
  }

  this->getLog().debug << "Calculating QM energy..." << Core::Log::endl;
  auto qmResults = Utils::CalculationRoutines::calculateWithCatch(*qmCalculator_, getLog(), "QM calculation failed");
  double qmEnergy = qmResults.get<Utils::Property::Energy>();
  this->getLog().debug << "Calculating MM energy..." << Core::Log::endl;
  Utils::Results mmResults;
  double mmEnergy;
  Utils::GradientCollection mmGradients;

  if (!calculateReducedQmMmEnergy_) {
    mmResults = Utils::CalculationRoutines::calculateWithCatch(*mmCalculator_, getLog(), "MM calculation failed");
    mmEnergy = mmResults.get<Utils::Property::Energy>();
    if (requiredProperties_.containsSubSet(Utils::Property::Gradients)) {
      mmGradients = mmResults.get<Utils::Property::Gradients>();
    }
  }
  else {
    prepareTermsForMmCalculator(true);
    mmResults = Utils::CalculationRoutines::calculateWithCatch(*mmCalculator_, getLog(), "MM calculation failed");
    mmEnergy = mmResults.get<Utils::Property::Energy>();
    if (requiredProperties_.containsSubSet(Utils::Property::Gradients)) {
      mmGradients = mmResults.get<Utils::Property::Gradients>();
    }
    prepareTermsForMmCalculator(false);
  }
  // Output:
  this->getLog().debug << std::setprecision(10) << "QM energy (Hartree): " << qmEnergy << Core::Log::endl;
  this->getLog().debug << "MM energy (Hartree): " << mmEnergy << Core::Log::endl;

  // Total energy:
  double totalEnergy = qmEnergy + mmEnergy;

  // Assemble results
  results_.set<Utils::Property::Description>(std::move(description));
  if (requiredProperties_.containsSubSet(Utils::Property::Energy)) {
    results_.set<Utils::Property::Energy>(totalEnergy);
  }
  if (requiredProperties_.containsSubSet(Utils::Property::Gradients)) {
    Utils::GradientCollection pcGradientsContributions;
    if (electrostaticEmbedding_ && mmAtomsLeft_)
      pcGradientsContributions = qmResults.get<Utils::Property::PointChargesGradients>();

    QmmmGradientsEvaluator gradEvaluator(qmResults.get<Utils::Property::Gradients>(), mmGradients,
                                         pcGradientsContributions, listOfQmAtoms_, mmBoundaryAtoms_,
                                         mmCalculator_->listsOfNeighbors(), structure_, qmRegion_);
    results_.set<Utils::Property::Gradients>(gradEvaluator.calculateQmmmGradients());
  }

  if (requiredProperties_.containsSubSet(Utils::Property::PartialHessian)) {
    QmmmHessianEvaluator hessianEvaluator(qmCalculator_, listOfQmAtoms_);
    results_.set<Utils::Property::PartialHessian>(hessianEvaluator.calculatePartialHessian());
  }
  if (requiredProperties_.containsSubSet(Utils::Property::AtomicCharges)) {
    auto totalCharges = mmResults.get<Utils::Property::AtomicCharges>();
    auto qmCharges = qmResults.get<Utils::Property::AtomicCharges>();
    int indexInQm = 0;
    for (const auto& indexInTotal : listOfQmAtoms_) {
      totalCharges[indexInTotal] = qmCharges[indexInQm++];
    }
    results_.set<Utils::Property::AtomicCharges>(totalCharges);
  }
  if (requiredProperties_.containsSubSet(Utils::Property::BondOrderMatrix)) {
    auto totalBondOrders = mmResults.get<Utils::Property::BondOrderMatrix>();
    const auto qmBondOrders = qmResults.get<Utils::Property::BondOrderMatrix>();
    int nQm = static_cast<int>(listOfQmAtoms_.size());
    for (int i = 0; i < nQm - 1; ++i) {
      auto totalIndexI = listOfQmAtoms_[i];
      for (int j = i + 1; j < nQm; ++j) {
        auto totalIndexJ = listOfQmAtoms_[j];
        totalBondOrders.setOrder(totalIndexI, totalIndexJ, qmBondOrders.getOrder(i, j));
      }
    }
    results_.set<Utils::Property::BondOrderMatrix>(totalBondOrders);
  }
  if (requiredProperties_.containsSubSet(Utils::Property::PartialEnergies)) {
    std::unordered_map<std::string, double> partialEnergies;
    partialEnergies.insert(std::make_pair("qm_energy", qmEnergy));
    partialEnergies.insert(std::make_pair("mm_energy", mmEnergy));
    if (!mmResults.has<Utils::Property::PartialEnergies>()) {
      throw std::runtime_error("MM calculator does not provide partial energies!");
    }
    for (const auto& partialResult : mmResults.get<Utils::Property::PartialEnergies>()) {
      partialEnergies.insert(partialResult);
    }
    if (qmResults.has<Utils::Property::PartialEnergies>()) {
      for (const auto& partialResult : qmResults.get<Utils::Property::PartialEnergies>()) {
        partialEnergies.insert(partialResult);
      }
    }
    results_.set<Utils::Property::PartialEnergies>(partialEnergies);
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

  handleElectrostaticEmbedding();
}

const Utils::PositionCollection& QmmmCalculator::getPositions() const {
  return structure_.getPositions();
}

void QmmmCalculator::setRequiredProperties(const Utils::PropertyList& requiredProperties) {
  // qmmm properties
  requiredProperties_ = requiredProperties;
  // qm properties
  auto qmProperties = requiredProperties_;
  // mm properties
  auto mmProperties = requiredProperties_;
  if (requiredProperties_.containsSubSet(Utils::Property::PartialHessian)) {
    mmProperties.removeProperty(Utils::Property::PartialHessian);
    qmProperties.removeProperty(Utils::Property::PartialHessian);
    qmProperties.addProperty(Utils::Property::Hessian);
  }
  mmCalculator_->setRequiredProperties(mmProperties);

  // This setting has to be recopied from the settings object, because it might have changed recently:
  electrostaticEmbedding_ = settings_->getBool(Utils::SettingsNames::electrostaticEmbedding);

  if (qmCalculator_) {
    if (requiredProperties_.containsSubSet(Utils::Property::PartialEnergies) &&
        !qmCalculator_->possibleProperties().containsSubSet(Utils::Property::PartialEnergies)) {
      qmProperties.removeProperty(Utils::Property::PartialEnergies);
    }
    if (electrostaticEmbedding_ && requiredProperties_.containsSubSet(Utils::Property::Gradients) && mmAtomsLeft_) {
      qmProperties.addProperty(Utils::Property::PointChargesGradients);
      qmCalculator_->setRequiredProperties(qmProperties);
    }
    else {
      qmCalculator_->setRequiredProperties(qmProperties);
    }
  }
}

Utils::PropertyList QmmmCalculator::getRequiredProperties() const {
  return requiredProperties_;
}

Utils::PropertyList QmmmCalculator::possibleProperties() const {
  Utils::PropertyList properties{Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::PartialHessian |
                                 Utils::Property::SuccessfulCalculation | Utils::Property::PartialEnergies};
  if (qmCalculator_ && qmCalculator_->possibleProperties().containsSubSet(Utils::Property::BondOrderMatrix)) {
    properties.addProperties(Utils::Property::BondOrderMatrix);
  }
  if (qmCalculator_ && qmCalculator_->possibleProperties().containsSubSet(Utils::Property::AtomicCharges)) {
    properties.addProperties(Utils::Property::AtomicCharges);
  }
  return properties;
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
      if (settings_->valueExists(key)) {
        mmCalculator_->settings().modifyValue(key, settings_->getValue(key));
      }
    }
  }

  // QM calculator
  if (qmCalculator_) {
    handleElectrostaticEmbedding();
    for (const auto& d : qmCalculator_->settings().getDescriptorCollection()) {
      if (!this->settings().valueExists(d.first))
        continue;
      if (d.first == Utils::ExternalQC::SettingsNames::pointChargesFile)
        continue;
      if (d.first == Utils::SettingsNames::mmCharges)
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
  auto log = getLog();
  if (this->mmCalculator_ != nullptr) {
    if (settings().getBool(SwooseUtilities::SettingsNames::silenceUnderlyingCalculators)) {
      Utils::CalculationRoutines::setLog(*mmCalculator_, log.error.operator bool(), log.warning.operator bool(), false);
    }
    else {
      this->mmCalculator_->setLog(log);
    }
  }
  if (this->qmCalculator_ != nullptr) {
    if (settings().getBool(SwooseUtilities::SettingsNames::silenceUnderlyingCalculators)) {
      Utils::CalculationRoutines::setLog(*qmCalculator_, log.error.operator bool(), log.warning.operator bool(), false);
    }
    else {
      this->qmCalculator_->setLog(log);
    }
  }
}

QmmmCalculator::~QmmmCalculator() {
  boost::filesystem::remove(pointChargesFilename_);
}

const Utils::Results& QmmmCalculator::ignoreQmCalculateImpl(std::string description) {
  this->getLog().debug << "\n\n";
  this->getLog().debug << "Ignoring QM calculation..." << Core::Log::nl << "Calculating MM energy..." << Core::Log::endl;
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

std::shared_ptr<MolecularMechanics::MolecularMechanicsCalculator> QmmmCalculator::getMolecularMechanicsCalculator() {
  return this->mmCalculator_;
}

std::shared_ptr<Core::Calculator> QmmmCalculator::getQuantumMechanicsCalculator() {
  return this->qmCalculator_;
}

} // namespace Qmmm
} // namespace Scine
