/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaffMolecularMechanicsCalculator.h"
#include "../Topology/IndexedStructuralTopologyCreator.h"
#include "GaffCalculatorSettings.h"
#include "GaffParameterDefaultsProvider.h"
#include "GaffParameterParser.h"
#include "GaffPotentialTermsGenerator.h"
#include <Core/Log.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondDetector.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <Utils/Math/FullSecondDerivativeCollection.h>

namespace Scine {
namespace MolecularMechanics {

std::string GaffMolecularMechanicsCalculator::name() const {
  return "GAFF";
}

GaffMolecularMechanicsCalculator::GaffMolecularMechanicsCalculator()
  : bondsEvaluator_(std::make_unique<BondsEvaluator>(structure_.getPositions())),
    anglesEvaluator_(std::make_unique<AnglesEvaluator>(structure_.getPositions())),
    dihedralsEvaluator_(std::make_unique<DihedralsEvaluator>(structure_.getPositions())),
    improperDihedralsEvaluator_(std::make_unique<DihedralsEvaluator>(structure_.getPositions())),
    electrostaticEvaluator_(std::make_unique<ElectrostaticEvaluator>(structure_.getPositions(), atomicCharges_)),
    lennardJonesEvaluator_(std::make_unique<LennardJonesEvaluator>(structure_.getPositions())) {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<GaffCalculatorSettings>();
  applySettings();
}

GaffMolecularMechanicsCalculator::GaffMolecularMechanicsCalculator(const GaffMolecularMechanicsCalculator& rhs)
  : CloneInterface(rhs), parameterFilePath_(rhs.parameterFilePath_), hessianMode_(rhs.hessianMode_) {
  this->bondsEvaluator_ = std::make_unique<BondsEvaluator>(structure_.getPositions());
  this->anglesEvaluator_ = std::make_unique<AnglesEvaluator>(structure_.getPositions());
  this->dihedralsEvaluator_ = std::make_unique<DihedralsEvaluator>(structure_.getPositions());
  this->improperDihedralsEvaluator_ = std::make_unique<DihedralsEvaluator>(structure_.getPositions());
  this->electrostaticEvaluator_ = std::make_unique<ElectrostaticEvaluator>(structure_.getPositions(), atomicCharges_);
  this->lennardJonesEvaluator_ = std::make_unique<LennardJonesEvaluator>(structure_.getPositions());

  this->requiredProperties_ = rhs.requiredProperties_;
  this->setLog(rhs.getLog());
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  applySettings();
  results_ = rhs.results_;
  this->listsOfNeighbors_ = rhs.listsOfNeighbors_; // has to be set before "setStructure"
  this->setParameters(rhs.parameters_);            // has to be set before "setStructure"
  structure_ = rhs.structure_;                     // has to be set last
  applySettings();
  if (structure_.size() != 0) {
    initialize();
  }
}

void GaffMolecularMechanicsCalculator::applySettings() {
  using namespace SwooseUtilities::SettingsNames;
  settings_->normalizeStringCases(); // convert all option names to lower case letters
  if (settings_->valid()) {
    applyCutoffDuringInitialization_ = settings_->getBool(applyCutoffDuringInitialization);
    nonCovalentCutoffRadius_ = settings_->getDouble(nonCovalentCutoffRadius) * Utils::Constants::bohr_per_angstrom;
    detectBondsWithCovalentRadii_ = settings_->getBool(detectBondsWithCovalentRadii);
    connectivityFilePath_ = settings_->getString(connectivityFilePath);
    onlyCalculateBondedContribution_ = settings_->getBool(onlyCalculateBondedContribution);
    printContributionsMolecularMechanics_ = settings_->getBool(printContributionsMolecularMechanics);
    atomicChargesFile_ = settings_->getString(gaffAtomicChargesFile);
    atomTypesFile_ = settings_->getString(gaffAtomTypesFile);

    std::string parameterFilePathFromSettings = settings_->getString(Utils::SettingsNames::parameterFilePath);
    if (parameterFilePath_ != parameterFilePathFromSettings) {
      parameterFilePath_ = parameterFilePathFromSettings;
      parameterFilePathHasBeenChanged_ = true;
    }
  }
  else {
    settings_->throwIncorrectSettings();
  }
}

void GaffMolecularMechanicsCalculator::setStructure(const Utils::AtomCollection& structure) {
  structure_ = structure;
  applySettings();
  initialize();
}

const Utils::Results& GaffMolecularMechanicsCalculator::calculate(std::string description) {
  applySettings();
  try {
    return calculateImpl(description);
  }
  catch (std::runtime_error& e) {
    throw Core::UnsuccessfulCalculationException(e.what());
  }
}

/*
 * Implementation of a calculation
 */
const Utils::Results& GaffMolecularMechanicsCalculator::calculateImpl(std::string description) {
  int nAtoms = structure_.size();
  double energy = 0.0;
  //  Initialize the atomic derivatives container
  Utils::AtomicSecondDerivativeCollection derivativesForBondedInteractions(nAtoms);
  derivativesForBondedInteractions.setZero();

  //  Initialize the full derivatives container
  int nAtomsInitialization = 0;
  if (!onlyCalculateBondedContribution_)
    nAtomsInitialization = nAtoms;
  const unsigned int derivativeOrder = (requiredProperties_.containsSubSet(Utils::Property::Hessian)) ? 2 : 1;

  Utils::DerivativeCollection fullDerivatives(nAtomsInitialization, derivativeOrder);
  fullDerivatives.setZero();
  double energyBonds = bondsEvaluator_->evaluate(derivativesForBondedInteractions);
  energy += energyBonds;
  double energyAngles = anglesEvaluator_->evaluate(derivativesForBondedInteractions);
  energy += energyAngles;
  double energyDihedrals = dihedralsEvaluator_->evaluate(derivativesForBondedInteractions);
  energy += energyDihedrals;
  double energyImproperDihedrals = improperDihedralsEvaluator_->evaluate(derivativesForBondedInteractions);
  energy += energyImproperDihedrals;
  double energyLennardJones = 0.0;
  double energyElectro = 0.0;
  if (!onlyCalculateBondedContribution_) {
    energyElectro = electrostaticEvaluator_->evaluate(fullDerivatives);
    energyLennardJones = lennardJonesEvaluator_->evaluate(fullDerivatives);
  }
  energy += energyLennardJones;
  energy += energyElectro;

  if (printContributionsMolecularMechanics_ && !hessianMode_) {
    this->getLog().output << Core::Log::nl << "Detailed output (unit: kcal/mol):" << Core::Log::endl;
    this->getLog().output << "-----------------------------------" << Core::Log::endl;
    this->getLog().output << "Bond energy: " << energyBonds * Utils::Constants::kCalPerMol_per_hartree << Core::Log::endl;
    this->getLog().output << "Angle energy: " << energyAngles * Utils::Constants::kCalPerMol_per_hartree << Core::Log::endl;
    this->getLog().output << "Dihedral energy: " << energyDihedrals * Utils::Constants::kCalPerMol_per_hartree
                          << Core::Log::endl;
    this->getLog().output << "Improper dihedral energy: " << energyImproperDihedrals * Utils::Constants::kCalPerMol_per_hartree
                          << Core::Log::endl;
    this->getLog().output << "Van der Waals (LJ) energy: " << energyLennardJones * Utils::Constants::kCalPerMol_per_hartree
                          << Core::Log::endl;
    this->getLog().output << "Electrostatic energy: " << energyElectro * Utils::Constants::kCalPerMol_per_hartree
                          << Core::Log::nl << Core::Log::endl;
  }

  // Assemble results
  results_.set<Utils::Property::Description>(std::move(description));
  if (requiredProperties_.containsSubSet(Utils::Property::Energy))
    results_.set<Utils::Property::Energy>(energy);
  if (requiredProperties_.containsSubSet(Utils::Property::Gradients)) {
    Utils::GradientCollection gradients(nAtoms, 3);
    // Add gradients from covalent contributions.
    for (int i = 0; i < nAtoms; ++i) {
      gradients.row(i) = Eigen::RowVector3d(derivativesForBondedInteractions[i].deriv());
    }
    // Add gradients from non-covalent contributions.
    if (!onlyCalculateBondedContribution_)
      gradients += fullDerivatives.getReferenceGradients();
    results_.set<Utils::Property::Gradients>(gradients);
  }
  if (requiredProperties_.containsSubSet(Utils::Property::Hessian)) {
    // Analytic Hessian for non-covalent contributions.
    const Utils::HessianMatrix& fullHessian = fullDerivatives.getHessianMatrix();
    // Add semi-numerical Hessian for covalent contributions.
    this->settings().modifyBool(SwooseUtilities::SettingsNames::onlyCalculateBondedContribution, true);
    this->hessianMode_ = true;
    Utils::NumericalHessianCalculator hessianCalculator(*this);
    Utils::Results numericalResult = hessianCalculator.calculate();
    this->hessianMode_ = false;

    if (!onlyCalculateBondedContribution_) {
      results_.set<Utils::Property::Hessian>(numericalResult.take<Utils::Property::Hessian>() + fullHessian);
      this->settings().modifyBool(SwooseUtilities::SettingsNames::onlyCalculateBondedContribution, false);
    }
    else {
      results_.set<Utils::Property::Hessian>(numericalResult.take<Utils::Property::Hessian>());
    }
  }
  if (requiredProperties_.containsSubSet(Utils::Property::AtomicCharges)) {
    results_.set<Utils::Property::AtomicCharges>(electrostaticEvaluator_->getAtomicCharges());
  }
  if (requiredProperties_.containsSubSet(Utils::Property::BondOrderMatrix)) {
    results_.set<Utils::Property::BondOrderMatrix>(
        SwooseUtilities::TopologyUtils::generateBondOrderMatrixFromListsOfNeighbors(listsOfNeighbors_));
  }
  if (requiredProperties_.containsSubSet(Utils::Property::PartialEnergies)) {
    std::unordered_map<std::string, double> partialEnergies;
    partialEnergies.insert(std::make_pair("bonds", energyBonds));
    partialEnergies.insert(std::make_pair("angles", energyAngles));
    partialEnergies.insert(std::make_pair("dihedral", energyDihedrals));
    partialEnergies.insert(std::make_pair("improper_dihedral", energyImproperDihedrals));
    partialEnergies.insert(std::make_pair("lennard_jones", energyLennardJones));
    partialEnergies.insert(std::make_pair("electrostatic", energyElectro));
    results_.set<Utils::Property::PartialEnergies>(partialEnergies);
  }
  results_.set<Utils::Property::SuccessfulCalculation>(true);

  return results_;
}

void GaffMolecularMechanicsCalculator::generatePotentialTerms(const std::string& parameterPath) {
  IndexedStructuralTopologyCreator topologyCreator(listsOfNeighbors_);
  auto topology = topologyCreator.calculateIndexedStructuralTopology();

  // Find out atom types
  GaffAtomTypeIdentifier atomTypeGenerator(structure_.getElements().size(), structure_.getElements(), listsOfNeighbors_,
                                           atomTypesFile_);
  auto atomTypes = atomTypeGenerator.getAtomTypes();
  if (!parametersHaveBeenSetInternally_ || parameterFilePathHasBeenChanged_) {
    if (parameterPath.empty()) {
      GaffParameterDefaultsProvider parameterProvider;
      this->getLog().output << "No parameter file was specified. Using the default GAFF parameters." << Core::Log::endl;
      parameters_ = *parameterProvider.getParameters();
    }
    else {
      GaffParameterParser parser(parameterPath);
      this->getLog().output << "Parsing the parameter file..." << Core::Log::endl;
      parameters_ = *parser.parseParameters();
      this->getLog().output << "Done." << Core::Log::nl << Core::Log::endl;
    }
    parametersHaveBeenSetInternally_ = true;
    parameterFilePathHasBeenChanged_ = false;
  }
  lennardJonesEvaluator_->setAtomTypesHolder(std::make_shared<AtomTypesHolder>(atomTypes));
  lennardJonesEvaluator_->setParameters(std::make_shared<GaffParameters>(parameters_));
  generatePotentialTerms(parameters_, topology, atomTypes);
}

void GaffMolecularMechanicsCalculator::generatePotentialTerms(const GaffParameters& parameters,
                                                              const IndexedStructuralTopology& topology,
                                                              const AtomTypesHolder& atomTypes) {
  // Get atomic charges
  if (atomicChargesFile_.empty()) {
    throw std::runtime_error("No atomic charges file provided! No electrostatic contributions will be evaluated. "
                             "Please provide an atomic charges file");
  }
  else {
    Eigen::MatrixXd chargesMatrix = Utils::csvToMatrix(atomicChargesFile_, ',');
    if (chargesMatrix.size() != structure_.size())
      throw std::runtime_error("The GAFF atomic charges file is invalid.");
    std::vector<double> charges(chargesMatrix.data(), chargesMatrix.data() + chargesMatrix.size());
    atomicCharges_ = std::move(charges);
  }

  // Generate the potential terms
  GaffPotentialTermsGenerator generator(structure_.size(), atomTypes, topology, parameters, structure_.getPositions(),
                                        this->getLog());

  bondsEvaluator_->setBondTerms(generator.getBondedTerms());
  anglesEvaluator_->setAngleTerms(generator.getAngleTerms());
  dihedralsEvaluator_->setDihedralTerms(generator.getDihedralTerms());
  improperDihedralsEvaluator_->setDihedralTerms(generator.getImproperDihedralTerms());
  if (!onlyCalculateBondedContribution_) {
    lennardJonesEvaluator_->setCutOffRadius(std::make_shared<double>(nonCovalentCutoffRadius_));
    lennardJonesEvaluator_->resetExclusions(structure_.size());
    lennardJonesEvaluator_->addExclusions(topology);
    lennardJonesEvaluator_->resetScaledInteractions(structure_.size());
    lennardJonesEvaluator_->addScaledInteractionPairs(topology);

    electrostaticEvaluator_->setCutOffRadius(std::make_shared<double>(nonCovalentCutoffRadius_));
    electrostaticEvaluator_->resetExclusions(structure_.size());
    electrostaticEvaluator_->addExclusions(topology);
    electrostaticEvaluator_->resetScaledInteractions(structure_.size());
    electrostaticEvaluator_->addScaledInteractionPairs(topology);
  }
}

void GaffMolecularMechanicsCalculator::initialize() {
  if (listsOfNeighbors_.empty()) {
    if (!detectBondsWithCovalentRadii_) {
      if (!connectivityFilePath_.empty()) {
        listsOfNeighbors_ = SwooseUtilities::ConnectivityFileHandler::readListsOfNeighbors(connectivityFilePath_);
        if (int(listsOfNeighbors_.size()) != structure_.size())
          throw std::runtime_error("The number of atoms in the provided connectivity file does not match the one of "
                                   "the molecular structure.");
      }
      else {
        throw std::runtime_error("The connectivity of the system could not be generated. Either specify a connectivity "
                                 "file or allow for the detection of bonds based on covalent radii (setting name: "
                                 "covalent_radii_bond_detection).");
      }
    }
    else {
      Utils::BondOrderCollection bondOrders = Utils::BondDetector::detectBonds(structure_);
      listsOfNeighbors_ =
          SwooseUtilities::TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(structure_.size(), bondOrders, 0.5);
    }
  }
  generatePotentialTerms(parameterFilePath_);
}

// If a Hessian is calculated, the calculator is cloned several times (if calculated in parallel),
// hence, parameters need to be stored as member.
void GaffMolecularMechanicsCalculator::setParameters(GaffParameters parameters) {
  parameters_ = std::move(parameters);
  parametersHaveBeenSetInternally_ = !parameters_.empty();
}

} // namespace MolecularMechanics
} // namespace Scine
