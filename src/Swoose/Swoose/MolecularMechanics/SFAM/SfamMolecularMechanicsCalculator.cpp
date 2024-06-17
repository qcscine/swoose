/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SfamMolecularMechanicsCalculator.h"
#include "../Parameters/DispersionC6Parameters.h"
#include "../Topology/IndexedStructuralTopologyCreator.h"
#include "SfamCalculatorSettings.h"
#include "SfamParameterParser.h"
#include "SfamPotentialTermsGenerator.h"
#include <Core/Log.h>
#include <Swoose/MMParametrization/ParametrizationUtils/ParameterFileWriter.h>
#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondDetector.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/Math/FullSecondDerivativeCollection.h>

namespace Scine {
namespace MolecularMechanics {

std::string SfamMolecularMechanicsCalculator::name() const {
  return "SFAM";
}

SfamMolecularMechanicsCalculator::SfamMolecularMechanicsCalculator()
  : bondsEvaluator_(std::make_unique<BondsEvaluator>(structure_.getPositions())),
    anglesEvaluator_(std::make_unique<AnglesEvaluator>(structure_.getPositions())),
    dihedralsEvaluator_(std::make_unique<DihedralsEvaluator>(structure_.getPositions())),
    improperDihedralsEvaluator_(std::make_unique<ImproperDihedralsEvaluator>(structure_.getPositions())),
    dispersionEvaluator_(std::make_unique<DispersionEvaluator>(structure_)),
    repulsionEvaluator_(std::make_unique<RepulsionEvaluator>(structure_)),
    electrostaticEvaluator_(std::make_unique<ElectrostaticEvaluator>(structure_.getPositions(), atomicCharges_)),
    hydrogenBondEvaluator_(std::make_unique<HydrogenBondEvaluator>(structure_, atomicCharges_)) {
  requiredProperties_ = Utils::Property::Energy;
  this->settings_ = std::make_unique<SfamCalculatorSettings>();
  applySettings();
}

SfamMolecularMechanicsCalculator::SfamMolecularMechanicsCalculator(const SfamMolecularMechanicsCalculator& rhs)
  : CloneInterface(rhs), parameterFilePath_(rhs.parameterFilePath_), hessianMode_(rhs.hessianMode_) {
  this->bondsEvaluator_ = std::make_unique<BondsEvaluator>(structure_.getPositions());
  this->anglesEvaluator_ = std::make_unique<AnglesEvaluator>(structure_.getPositions());
  this->dihedralsEvaluator_ = std::make_unique<DihedralsEvaluator>(structure_.getPositions());
  this->improperDihedralsEvaluator_ = std::make_unique<ImproperDihedralsEvaluator>(structure_.getPositions());
  this->dispersionEvaluator_ = std::make_unique<DispersionEvaluator>(structure_);
  this->repulsionEvaluator_ = std::make_unique<RepulsionEvaluator>(structure_);
  this->electrostaticEvaluator_ = std::make_unique<ElectrostaticEvaluator>(structure_.getPositions(), atomicCharges_);
  this->hydrogenBondEvaluator_ = std::make_unique<HydrogenBondEvaluator>(structure_, atomicCharges_);

  this->requiredProperties_ = rhs.requiredProperties_;
  this->setLog(rhs.getLog());
  auto valueCollection = dynamic_cast<const Utils::UniversalSettings::ValueCollection&>(rhs.settings());
  this->settings_ =
      std::make_unique<Utils::Settings>(Utils::Settings(valueCollection, rhs.settings().getDescriptorCollection()));
  applySettings();
  results_ = rhs.results_;
  this->setListsOfNeighbors(rhs.listsOfNeighbors_); // has to be set before "setStructure"
  this->setParameters(rhs.parameters_);             // has to be set before "setStructure"
  structure_ = rhs.structure_;                      // has to be set last
  applySettings();
  if (structure_.size() != 0) {
    initialize();
  }
}

void SfamMolecularMechanicsCalculator::applySettings() {
  using namespace SwooseUtilities::SettingsNames;
  settings_->normalizeStringCases(); // convert all option names to lower case letters
  if (settings_->valid()) {
    applyCutoffDuringInitialization_ = settings_->getBool(applyCutoffDuringInitialization);
    includeHydrogenBonds_ = settings_->getBool(hydrogenBondCorrection);
    nonCovalentCutoffRadius_ = settings_->getDouble(nonCovalentCutoffRadius) * Utils::Constants::bohr_per_angstrom;
    detectBondsWithCovalentRadii_ = settings_->getBool(detectBondsWithCovalentRadii);
    connectivityFilePath_ = settings_->getString(connectivityFilePath);
    onlyCalculateBondedContribution_ = settings_->getBool(onlyCalculateBondedContribution);
    printContributionsMolecularMechanics_ = settings_->getBool(printContributionsMolecularMechanics);
    std::string sfamAtomTypeLevelString = settings_->getString(sfamAtomTypeLevel);
    sfamAtomTypeLevel_ = SfamAtomTypeIdentifier::generateSfamAtomTypeLevelFromString(sfamAtomTypeLevelString);

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

void SfamMolecularMechanicsCalculator::setStructure(const Utils::AtomCollection& structure) {
  structure_ = structure;
  applySettings();
  initialize();
}

void SfamMolecularMechanicsCalculator::modifyPositions(Utils::PositionCollection newPositions) {
  structure_.setPositions(newPositions);
  IndexedStructuralTopologyCreator topologyCreator(listsOfNeighbors_);
  if (includeHydrogenBonds_) {
    topologyCreator.addHydrogenBondsToIndexedStructuralTopology(topology_, structure_);
    SfamPotentialTermsGenerator generator(structure_.size(), atomTypes_, topology_, parameters_,
                                          structure_.getPositions(), nonCovalentCutoffRadius_, this->getLog());
    hydrogenBondEvaluator_->setHydrogenBondTerms(generator.getHydrogenBondTerms());
  }
}

const Utils::Results& SfamMolecularMechanicsCalculator::calculate(std::string description) {
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
const Utils::Results& SfamMolecularMechanicsCalculator::calculateImpl(std::string description) {
  int nAtoms = structure_.size();
  double energy = 0.0;

  //  Initialize the atomic derivatives container
  Utils::AtomicSecondDerivativeCollection derivativesForBondedInteractions(nAtoms);
  derivativesForBondedInteractions.setZero();

  //  Initialize the full derivatives container
  int nAtomsInitialization = 0;
  if (!onlyCalculateBondedContribution_)
    nAtomsInitialization = nAtoms;
  Utils::FullSecondDerivativeCollection fullDerivatives(nAtomsInitialization);
  fullDerivatives.setZero();

  Utils::Dftd3::Dftd3 d3;
  auto d3Pointer = std::make_shared<Utils::Dftd3::Dftd3>(d3);
  Eigen::MatrixXd R0(nAtoms, nAtoms);

  double energyBonds = bondsEvaluator_->evaluate(derivativesForBondedInteractions);
  energy += energyBonds;
  double energyAngles = anglesEvaluator_->evaluate(derivativesForBondedInteractions);
  energy += energyAngles;
  double energyDihedrals = dihedralsEvaluator_->evaluate(derivativesForBondedInteractions);
  energy += energyDihedrals;
  double energyImproperDihedrals = improperDihedralsEvaluator_->evaluate(derivativesForBondedInteractions);
  energy += energyImproperDihedrals;

  double energyHydrogenBonds = 0.0;
  if (includeHydrogenBonds_) {
    energyHydrogenBonds = hydrogenBondEvaluator_->evaluate(derivativesForBondedInteractions);
  }
  energy += energyHydrogenBonds;

  double energyDispersion = 0.0;
  double energyRepulsion = 0.0;
  double energyElectro = 0.0;
  if (!onlyCalculateBondedContribution_) {
    energyElectro = electrostaticEvaluator_->evaluate(fullDerivatives);
    energyDispersion = dispersionEvaluator_->evaluate(fullDerivatives, d3Pointer, R0);
    energyRepulsion = repulsionEvaluator_->evaluate(fullDerivatives, R0);
  }
  energy += energyDispersion;
  energy += energyRepulsion;
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
    this->getLog().output << "Hydrogen Bond energy: " << energyHydrogenBonds * Utils::Constants::kCalPerMol_per_hartree
                          << Core::Log::endl;
    this->getLog().output << "Dispersion energy: " << energyDispersion * Utils::Constants::kCalPerMol_per_hartree
                          << Core::Log::endl;
    this->getLog().output << "Repulsion energy: " << energyRepulsion * Utils::Constants::kCalPerMol_per_hartree
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
    Utils::Results numericalResult;
    if (atomsToConsiderForHessian_.empty())
      numericalResult = hessianCalculator.calculate();
    else // In this case the Hessian is not the true Hessian of the system (be careful)
      numericalResult = hessianCalculator.calculate(atomsToConsiderForHessian_);
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
    partialEnergies.insert(std::make_pair("h_bonds", energyHydrogenBonds));
    partialEnergies.insert(std::make_pair("dispersion", energyDispersion));
    partialEnergies.insert(std::make_pair("repulsion", energyRepulsion));
    partialEnergies.insert(std::make_pair("electrostatic", energyElectro));
    results_.set<Utils::Property::PartialEnergies>(partialEnergies);
  }
  results_.set<Utils::Property::SuccessfulCalculation>(true);

  return results_;
}

void SfamMolecularMechanicsCalculator::generatePotentialTerms(const std::string& parameterPath) {
  IndexedStructuralTopologyCreator topologyCreator(listsOfNeighbors_);
  topology_ = topologyCreator.calculateIndexedStructuralTopology();
  if (includeHydrogenBonds_)
    topologyCreator.addHydrogenBondsToIndexedStructuralTopology(topology_, structure_);

  // Find out atom types
  SfamAtomTypeIdentifier atomTypeGenerator(structure_.getElements().size(), structure_.getElements(), listsOfNeighbors_);
  atomTypes_ = atomTypeGenerator.getAtomTypes(sfamAtomTypeLevel_);

  if (!parametersHaveBeenSetInternally_ || parameterFilePathHasBeenChanged_) {
    if (parameterPath.empty())
      throw std::runtime_error("No parameter file was specified.");

    SfamParameterParser parser(parameterPath, atomTypes_);
    this->getLog().output << "Parsing the parameter file..." << Core::Log::endl;
    parameters_ = *parser.parseParameters();
    parametersHaveBeenSetInternally_ = true;
    parameterFilePathHasBeenChanged_ = false;
    this->getLog().output << "Done." << Core::Log::nl << Core::Log::endl;
    // Sanity check
    if (!parameters_.sanityCheck(atomTypes_)) {
      // Try assigning new C6 coefficients and check again:
      DispersionC6Parameters::fillC6MatrixForCurrentStructure(parameters_, structure_, atomTypes_);
      if (!parameters_.sanityCheck(atomTypes_))
        throw std::runtime_error("The parameters, which were read from a file, are not valid!");

      // If check successful, then update parameters file with C6 parameters
      MMParametrization::ParameterFileWriter::writeSfamParametersToFile(parameterPath, parameters_);
      this->getLog().output << "As no C6 coefficients were provided, they were calculated on the fly and the "
                               "parameter file was updated accordingly.\n"
                            << Core::Log::endl;
    }
  }

  generatePotentialTerms(parameters_, topology_, atomTypes_);
}

void SfamMolecularMechanicsCalculator::generatePotentialTerms(const SfamParameters& parameters,
                                                              const IndexedStructuralTopology& topology,
                                                              const AtomTypesHolder& atomTypes) {
  // Get charges from parameters
  atomicCharges_ = parameters.getChargesForEachAtom(atomTypes);
  // Apply non-covalent global parameters
  auto nonCovalentParameters = parameters.getNonCovalentParameters();
  electrostaticEvaluator_->setScalingFactor(nonCovalentParameters.back());
  nonCovalentParameters.pop_back();
  repulsionEvaluator_->setBetaRepulsion(nonCovalentParameters.back());
  nonCovalentParameters.pop_back();
  dispersionEvaluator_->setD3Parameters(nonCovalentParameters);

  // Generate the potential terms
  SfamPotentialTermsGenerator generator(structure_.size(), atomTypes, topology, parameters, structure_.getPositions(),
                                        nonCovalentCutoffRadius_, this->getLog());

  bondsEvaluator_->setBondTerms(generator.getBondedTerms());
  anglesEvaluator_->setAngleTerms(generator.getAngleTerms());
  dihedralsEvaluator_->setDihedralTerms(generator.getDihedralTerms());
  improperDihedralsEvaluator_->setImproperDihedralTerms(generator.getImproperDihedralTerms());
  hydrogenBondEvaluator_->setHydrogenBondTerms(generator.getHydrogenBondTerms());
  if (!onlyCalculateBondedContribution_) {
    dispersionEvaluator_->setDispersionTerms(generator.getDispersionTerms(applyCutoffDuringInitialization_));
    repulsionEvaluator_->setRepulsionTerms(generator.getRepulsionTerms(applyCutoffDuringInitialization_));
    electrostaticEvaluator_->setElectrostaticTerms(generator.getElectrostaticTerms(applyCutoffDuringInitialization_));
  }
}

void SfamMolecularMechanicsCalculator::initialize() {
  for (const auto& element : structure_.getElements()) {
    if (Utils::ElementInfo::Z(element) > 94)
      throw std::runtime_error("This MM method is only defined for elements up to element 94.");
  }

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

// This function is only used by the friend class UpdateFunctionManager in the MM parametrization algorithm.
void SfamMolecularMechanicsCalculator::setListsOfNeighbors(std::vector<std::list<int>> listsOfNeighbors) {
  listsOfNeighbors_ = std::move(listsOfNeighbors);
}

// If a Hessian is calculated, the calculator is cloned several times (if calculated in parallel),
// hence, parameters need to be stored as member.
void SfamMolecularMechanicsCalculator::setParameters(SfamParameters parameters) {
  // Only set the parameters if there is actually something in it
  if (!parameters.getNonCovalentParameters().empty()) {
    parameters_ = std::move(parameters);
    parametersHaveBeenSetInternally_ = true;
  }
}

void SfamMolecularMechanicsCalculator::setAtomsToConsiderForHessian(std::vector<int> atomsToConsiderForHessian) {
  atomsToConsiderForHessian_ = std::move(atomsToConsiderForHessian);
}

} // namespace MolecularMechanics
} // namespace Scine
