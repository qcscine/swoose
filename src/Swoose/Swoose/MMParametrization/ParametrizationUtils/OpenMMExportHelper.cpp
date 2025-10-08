/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OpenMMExportHelper.h"
#include "Swoose/Utilities/SettingsNames.h"
#include <Swoose/MolecularMechanics/Interactions/HydrogenBondParameters.h>
#include <Swoose/MolecularMechanics/Interactions/RepulsionParameters.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/ChemicalFileFormats/PdbStreamHandler.h>
#include <Utils/Settings.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace MMParametrization {

OpenMMExportHelper::OpenMMExportHelper(std::shared_ptr<Utils::Settings> settings, ParametrizationData& data)
  : settings_(settings),
    data_(data),
    parameters_(data.parameters),
    nAtoms_(data_.fullStructure.size()),
    xmlFileName_(settings_->getString(SwooseUtilities::SettingsNames::sfamOpenMMFileName)) {
  initializeUniqueAtomTypes();
}

/*
 * At the moment it does not look like that we actually need to export the topology.
 */
void OpenMMExportHelper::exportTopology() {
  Utils::PdbStreamHandler pdbFileHandler;
  std::ofstream pdbFile("sfam.pdb");
  std::map<std::string, std::string> atomTypeToResidueName;
  for (unsigned int i = 0; i < uniqueAtomTypes_.size(); ++i) {
    atomTypeToResidueName.insert({uniqueAtomTypes_[i], "R" + std::to_string(i)});
  }
  Utils::ResidueCollection residues;
  for (unsigned int i = 0; i < nAtoms_; ++i) {
    const auto atomType = data_.atomTypes.getAtomType(int(i));
    residues.emplace_back(Utils::ResidueInformation({atomTypeToResidueName[atomType], "", "A", i + 1}));
  }
  data_.fullStructure.setResidues(residues);
  pdbFileHandler.write(pdbFile, "pdb", data_.fullStructure, data_.bondOrders,
                       "Structure and topology used during SFAM parametrization.");
}

void OpenMMExportHelper::exportSfamForceField() {
  std::cout << " == Exporting Force Field as XML file " << xmlFileName_ << " == " << std::endl;
  std::ofstream xmlFile(xmlFileName_);

  xmlFile << "<ForceField>\n";
  writeHeader(xmlFile);
  writeHarmonicBondTerms(xmlFile);
  writeAngleTerms(xmlFile);
  writeDihedralTerms(xmlFile);
  writeImproperDihedralTerms(xmlFile);
  writeNonbondedTerms(xmlFile);
  writeCoulombTerm(xmlFile);
  writeHydrogenBondTerms(xmlFile);
  xmlFile << "</ForceField>\n";
}

void OpenMMExportHelper::writeHeader(std::ofstream& xmlFile) {
  writeAtomTypes(xmlFile);
  writeResidues(xmlFile);
}

void OpenMMExportHelper::initializeUniqueAtomTypes() {
  uniqueAtomTypes_.clear();
  MolecularMechanics::RepulsionParameterHelper helper;
  for (unsigned int i = 0; i < nAtoms_; i++) {
    auto atomType = data_.atomTypes.getAtomType(int(i));
    auto elementType = data_.fullStructure.at(int(i)).getElementType();

    if (std::find(uniqueAtomTypes_.begin(), uniqueAtomTypes_.end(), atomType) == uniqueAtomTypes_.end()) {
      uniqueAtomTypes_.push_back(atomType);
      effectiveChargesHolder_.push_back(helper.getEffectiveCharge(elementType));
      elementTypes_.push_back(elementType);
    }
  }
}

void OpenMMExportHelper::writeAtomTypes(std::ofstream& xmlFile) {
  xmlFile << " <AtomTypes>\n";
  for (unsigned int i = 0; i < uniqueAtomTypes_.size(); ++i) {
    const std::string atomType = uniqueAtomTypes_[i];
    const auto element = Utils::ElementInfo::symbol(elementTypes_[i]);
    const auto mass = Utils::ElementInfo::mass(elementTypes_[i]);
    xmlFile << "  <Type name=" << '"' << atomType << '"';
    xmlFile << " class=" << '"' << element << '"';
    xmlFile << " element=" << '"' << element << '"';
    xmlFile << " mass=" << '"' << std::setprecision(8) << mass << '"' << "/>\n";
  }
  xmlFile << " </AtomTypes>\n";
}

void OpenMMExportHelper::writeResidues(std::ofstream& xmlFile) {
  // SFAM does not use any residues. Therefore, we have one residue per atom type.
  xmlFile << " <Residues>\n";
  const auto& charges = parameters_.getCharges();
  for (unsigned int i = 0; i < uniqueAtomTypes_.size(); ++i) {
    const std::string& atomType = uniqueAtomTypes_[i];
    const std::string& element = Utils::ElementInfo::symbol(elementTypes_[i]);
    const std::string atomName = element;
    xmlFile << "  <Residue name=" << '"' << "R" << i << '"' << ">\n";
    xmlFile << "    <Atom name=" << '"' << atomName << '"';
    xmlFile << " type=" << '"' << atomType << '"';
    xmlFile << " charge=" << '"' << charges.at(atomType) << '"' << "/>\n";
    xmlFile << "    <ExternalBond atomName=" << '"' << atomName << '"' << "/>\n";
    xmlFile << "  </Residue\n>";
  }
  xmlFile << " </Residues>\n";
}

void OpenMMExportHelper::writeHarmonicBondTerms(std::ofstream& xmlFile) {
  if (!parameters_.getBonds().empty()) {
    xmlFile << " <HarmonicBondForce>\n";
    for (const auto& bond : parameters_.getBonds()) {
      xmlFile << "  <Bond";
      xmlFile << " type1=" << '"' << bond.first.a1 << '"';
      xmlFile << " type2=" << '"' << bond.first.a2 << '"';
      // The factor of 0.1 is nanometer_to_angstrom conversion
      xmlFile << " length=" << '"' << bond.second.getEquilibriumBondLength() * 0.1 << '"';
      xmlFile << " k=" << '"' << bond.second.getForceConstant() * Utils::Constants::joule_per_calorie * 100.0 << '"' << "/>\n";
    }
    xmlFile << " </HarmonicBondForce>\n";
  }
}

void OpenMMExportHelper::writeAngleTerms(std::ofstream& xmlFile) {
  if (!parameters_.getAngles().empty()) {
    xmlFile << " <HarmonicAngleForce>\n";
    for (const auto& angle : parameters_.getAngles()) {
      xmlFile << "  <Angle";
      xmlFile << " type1=" << '"' << angle.first.a1 << '"';
      xmlFile << " type2=" << '"' << angle.first.a2 << '"';
      xmlFile << " type3=" << '"' << angle.first.a3 << '"';
      xmlFile << " angle=" << '"' << angle.second.getEquilibriumAngle() * Utils::Constants::rad_per_degree << '"';
      xmlFile << " k=" << '"' << angle.second.getForceConstant() * Utils::Constants::joule_per_calorie << '"' << "/>\n";
    }
    xmlFile << " </HarmonicAngleForce>\n";
  }
}

void OpenMMExportHelper::writeDihedralTerms(std::ofstream& xmlFile) {
  if (!parameters_.getDihedrals().empty()) {
    xmlFile << " <CustomTorsionForce energy=";
    const std::string formula = "halfBarrierHeight*(1-cos(n*theta-phase))";
    xmlFile << '"' << formula << '"';
    xmlFile << ">\n";
    xmlFile << "  <PerTorsionParameter name=" << '"' << "halfBarrierHeight" << '"' << "/>\n";
    xmlFile << "  <PerTorsionParameter name=" << '"' << "n" << '"' << "/>\n";
    xmlFile << "  <PerTorsionParameter name=" << '"' << "phase" << '"' << "/>\n";
    for (const auto& dihedral : parameters_.getDihedrals()) {
      xmlFile << "  <Proper ";
      const auto type1 = (dihedral.first.a1 != "X") ? dihedral.first.a1 : "";
      const auto type2 = (dihedral.first.a2 != "X") ? dihedral.first.a2 : "";
      const auto type3 = (dihedral.first.a3 != "X") ? dihedral.first.a3 : "";
      const auto type4 = (dihedral.first.a4 != "X") ? dihedral.first.a4 : "";
      xmlFile << " type1=" << '"' << type1 << '"';
      xmlFile << " type2=" << '"' << type2 << '"';
      xmlFile << " type3=" << '"' << type3 << '"';
      xmlFile << " type4=" << '"' << type4 << '"';
      xmlFile << " n=" << '"' << dihedral.second.getPeriodicity() << '"';
      xmlFile << " phase1=" << '"' << dihedral.second.getPhaseShift() * Utils::Constants::rad_per_degree << '"';
      xmlFile << " halfBarrierHeight=" << '"'
              << dihedral.second.getHalfBarrierHeight() * Utils::Constants::joule_per_calorie << '"' << "/>\n";
    }
    xmlFile << " </CustomTorsionForce>\n";
  }
}
void OpenMMExportHelper::writeImproperDihedralTerms(std::ofstream& xmlFile) {
  if (!parameters_.getImproperDihedrals().empty()) {
    xmlFile << " <CustomTorsionForce energy=";
    const std::string formula = "k*(cos(theta)-cos(phase))";
    xmlFile << '"' << formula << '"';
    xmlFile << ">\n";
    xmlFile << "  <PerTorsionParameter name=" << '"' << "k" << '"' << "/>\n";
    xmlFile << "  <PerTorsionParameter name=" << '"' << "phase" << '"' << "/>\n";
    for (const auto& improperDihedral : parameters_.getImproperDihedrals()) {
      xmlFile << "  <Improper ";
      const auto type1 = (improperDihedral.first.ac != "X") ? improperDihedral.first.ac : "";
      const auto type2 = (improperDihedral.first.a2 != "X") ? improperDihedral.first.a2 : "";
      const auto type3 = (improperDihedral.first.a3 != "X") ? improperDihedral.first.a3 : "";
      const auto type4 = (improperDihedral.first.a4 != "X") ? improperDihedral.first.a4 : "";
      xmlFile << " type1=" << '"' << type1 << '"';
      xmlFile << " type2=" << '"' << type2 << '"';
      xmlFile << " type3=" << '"' << type3 << '"';
      xmlFile << " type4=" << '"' << type4 << '"';
      xmlFile << " k=" << '"' << improperDihedral.second.getForceConstant() * Utils::Constants::joule_per_calorie << '"';
      xmlFile << " phase=" << '"' << improperDihedral.second.getEquilibriumAngle() * Utils::Constants::rad_per_degree
              << '"' << "/>\n";
    }
    xmlFile << " </CustomTorsionForce>\n";
  }
}

void OpenMMExportHelper::writeCoulombTerm(std::ofstream& xmlFile) {
  xmlFile << " <NonbondedForce bondCutoff=" << '"' << 3 << '"';
  xmlFile << " coulomb14scale=" << '"' << 1.0 << '"';
  xmlFile << " lj14scale=" << '"' << 0.0 << '"' << ">\n";

  const auto& charges = parameters_.getCharges();
  for (const auto& atomType : uniqueAtomTypes_) {
    xmlFile << "  <Atom type=" << '"' << atomType << '"';
    xmlFile << " charge=" << '"' << charges.at(atomType) << '"';
    xmlFile << " epsilon=" << '"' << 0.0 << '"';
    xmlFile << " sigma=" << '"' << 1.0 << '"' << "/>\n";
    ;
  }
  xmlFile << " </NonbondedForce>\n";
}

void OpenMMExportHelper::writeNonbondedTerms(std::ofstream& xmlFile) {
  // define the functional form
  xmlFile << " <CustomNonbondedForce bondCutoff=" << '"' << 3 << '"' << "\n";
  xmlFile << "  energy=" << '"' << "scale*(-dispTerm1-dispTerm2+repulsionTerm);\n";
  xmlFile << "   dispTerm1=C6(i1, i2)/(dist^6 + fdamp^6);\n";
  xmlFile << "   dispTerm2=Cs8(i1, i2)/(dist^8 + fdamp^8);\n";
  xmlFile << "   fdamp=a1*R0(i1, i2)+a2;\n";
  xmlFile << "   repulsionTerm=Zeff1*Zeff2/dist * exp(-beta*dist/R0(i1, i2));\n";
  xmlFile << "   dist = r*bohr_per_nm;" << '"' << ">\n";
  xmlFile << "  <GlobalParameter name=" << '"' << "a1" << '"';
  xmlFile << " defaultValue=" << '"' << parameters_.getNonCovalentParameters().at(0) << '"' << "/>\n";
  xmlFile << "  <GlobalParameter name=" << '"' << "a2" << '"';
  xmlFile << " defaultValue=" << '"' << parameters_.getNonCovalentParameters().at(2) << '"' << "/>\n";
  xmlFile << "  <GlobalParameter name=" << '"' << "beta" << '"';
  xmlFile << " defaultValue=" << '"' << parameters_.getNonCovalentParameters().at(3) << '"' << "/>\n";
  xmlFile << "  <GlobalParameter name=" << '"' << "scale" << '"';
  xmlFile << " defaultValue=" << '"' << Utils::Constants::kJPerMol_per_hartree << '"' << "/>\n";
  xmlFile << "  <GlobalParameter name=" << '"' << "bohr_per_nm" << '"';
  xmlFile << " defaultValue=" << '"' << Utils::Constants::bohr_per_angstrom * 10.0 << '"' << "/>\n";
  xmlFile << "  <PerParticleParameter name=" << '"' << "Zeff" << '"' << "/>\n";
  xmlFile << "  <PerParticleParameter name=" << '"' << "i" << '"' << "/>\n";

  xmlFile << "  <Function name=" << '"' << "C6" << '"';
  xmlFile << " type=" << '"' << "Discrete2D" << '"';
  xmlFile << " xsize=" << '"' << uniqueAtomTypes_.size() << '"';
  xmlFile << " ysize=" << '"' << uniqueAtomTypes_.size() << '"' << ">\n";
  xmlFile << parameters_.getC6Matrix() << "\n";
  xmlFile << "  </Function>\n";

  xmlFile << "  <Function name=" << '"' << "Cs8" << '"';
  xmlFile << " type=" << '"' << "Discrete2D" << '"';
  xmlFile << " xsize=" << '"' << uniqueAtomTypes_.size() << '"';
  xmlFile << " ysize=" << '"' << uniqueAtomTypes_.size() << '"' << ">\n";
  xmlFile << parameters_.getC8Matrix() * parameters_.getNonCovalentParameters().at(1) << "\n";
  xmlFile << "  </Function>\n";

  // Calculate R0 from C6 and C8
  const Eigen::MatrixXf R0 = calculateR0Matrix();
  xmlFile << "  <Function name=" << '"' << "R0" << '"';
  xmlFile << " type=" << '"' << "Discrete2D" << '"';
  xmlFile << " xsize=" << '"' << uniqueAtomTypes_.size() << '"';
  xmlFile << " ysize=" << '"' << uniqueAtomTypes_.size() << '"' << ">\n";
  xmlFile << R0 << "\n";
  xmlFile << "  </Function>\n";

  const auto indexMap = parameters_.getC6IndicesMap();
  for (unsigned int i = 0; i < uniqueAtomTypes_.size(); i++) {
    const std::string atomType = uniqueAtomTypes_.at(i);
    xmlFile << "  <Atom type=" << '"' << atomType << '"';
    xmlFile << " i=" << '"' << indexMap.at(atomType) << '"';
    xmlFile << " Zeff=" << '"' << effectiveChargesHolder_.at(i) << '"' << "/>\n";
  }

  xmlFile << " </CustomNonbondedForce>\n";
}

void OpenMMExportHelper::writeHydrogenBondTerms(std::ofstream& xmlFile) {
  const double kappa1 = 10.0;
  const double kappa2 = 5.0;
  const double nanometerPerBohr = Utils::Constants::angstrom_per_bohr / 10.0;
  const double nanometerPerBohrCubed = nanometerPerBohr * nanometerPerBohr * nanometerPerBohr;
  xmlFile << " <CustomHbondForce";

  xmlFile << " particlesPerDonor=" << '"' << "2" << '"';    // D-H
  xmlFile << " particlesPerAcceptor=" << '"' << "1" << '"'; // A
  xmlFile << " bondCutoff=" << '"' << "3" << '"';

  xmlFile << "  energy=" << '"' << "-fdamp*(combA+combD)/(r^3);\n";
  xmlFile << "   fdamp=firstTerm*secondTerm;\n";
  xmlFile << "   firstTerm=1.0/(1.0+(r/tilde_r)^(12));\n";
  xmlFile << "   secondTerm=(0.5*(cos(phi)+1))^6;\n";
  xmlFile << "   r = distance(a1, d2);\n";
  xmlFile << "   phi = angle(d2, d1, a1) + pi;" << '"' << ">\n";

  xmlFile << "  <GlobalParameter name=" << '"' << "tilde_r" << '"' << " defaultValue=" << '"' << 0.4 << '"' << "/>\n"; // in nm
  xmlFile << "  <GlobalParameter name=" << '"' << "pi" << '"' << " defaultValue=" << '"' << Utils::Constants::pi << '"'
          << "/>\n"; // in nm

  xmlFile << "  <PerDonorParameter name=" << '"' << "combD" << '"' << "/>\n";
  xmlFile << "  <PerAcceptorParameter name=" << '"' << "combA" << '"' << "/>\n";

  detectDonorAndAcceptorGroups();

  for (auto& donor : donors_) {
    const double expKD = exp(-kappa1 * donor.chargeD);
    // The atom ordering is important here. The first atom must be the H-atom. Otherwise, the H-bond will be inverted.
    xmlFile << "  <Donor";
    xmlFile << " type1=" << '"' << donor.hAtomType << '"';
    xmlFile << " type2=" << '"' << donor.dAtomType << '"';
    xmlFile << " combD=" << '"'
            << donor.kHBD * expKD / (expKD + kappa2) * Utils::Constants::kJPerMol_per_hartree * nanometerPerBohrCubed
            << '"'; // in kJ/mol * nm^3
    xmlFile << "/>\n";
  }

  for (auto& acceptor : acceptors_) {
    const double expKA = exp(-kappa1 * acceptor.chargeA);
    xmlFile << "  <Acceptor";
    xmlFile << " type1=" << '"' << acceptor.aAtomType << '"';
    xmlFile << " combA=" << '"'
            << acceptor.kHBA * expKA / (expKA + kappa2) * Utils::Constants::kJPerMol_per_hartree * nanometerPerBohrCubed
            << '"';
    xmlFile << "/>\n";
  }

  xmlFile << " </CustomHbondForce>\n";
}

Eigen::MatrixXf OpenMMExportHelper::calculateR0Matrix() {
  Eigen::MatrixXf R0 = Eigen::MatrixXf::Zero(uniqueAtomTypes_.size(), uniqueAtomTypes_.size());
  const auto& C8 = parameters_.getC8Matrix();
  const auto& C6 = parameters_.getC6Matrix();
  for (unsigned int i = 0; i < uniqueAtomTypes_.size(); i++) {
    for (unsigned int j = 0; j < uniqueAtomTypes_.size(); j++) {
      const double& c6ij = C6(i, j);
      R0(i, j) = (c6ij > 1e-12) ? std::sqrt(C8(i, j) / C6(i, j)) : 1.0;
    }
  }
  return R0;
}

std::pair<bool, int> OpenMMExportHelper::isBondedToHydrogenAtom(int atomIndex) {
  auto listOfNeighbors = data_.listsOfNeighbors.at(atomIndex);
  for (auto& neighbor : listOfNeighbors) {
    if (data_.fullStructure.at(neighbor).getElementType() == Utils::ElementType::H) {
      return std::make_pair(true, neighbor);
    }
  }
  return std::make_pair(false, -1);
}

void OpenMMExportHelper::detectDonorAndAcceptorGroups() {
  using namespace MolecularMechanics::HydrogenBondHelper;
  for (int i = 0; i < int(nAtoms_); i++) {
    auto atomType = data_.atomTypes.getAtomType(i);
    auto elementType = data_.fullStructure.at(i).getElementType();
    const double charge = parameters_.getCharges().find(atomType)->second;
    auto hbTypes = MolecularMechanics::HydrogenBondHelper::vectorOfDonorOrAcceptorElements_;

    if (std::find(hbTypes.begin(), hbTypes.end(), elementType) != hbTypes.end()) {
      const double interactionStrength =
          MolecularMechanics::HydrogenBondParameters::getInteractionStrengthConstants(elementType);
      if (isBondedToHydrogenAtom(i).first) {
        const int correspondingHIndex = isBondedToHydrogenAtom(i).second;
        const Donor donor(atomType, data_.atomTypes.getAtomType(correspondingHIndex), interactionStrength, charge);
        auto it = std::find_if(donors_.begin(), donors_.end(), [donor](Donor d) { return donor == d; });
        if (it == donors_.end()) {
          donors_.push_back(donor);
        }
      }
      const Acceptor acceptor(atomType, interactionStrength, charge);
      auto it = std::find_if(acceptors_.begin(), acceptors_.end(), [acceptor](Acceptor a) { return acceptor == a; });
      if (it == acceptors_.end()) {
        acceptors_.push_back(acceptor);
      }
    }
  }
}

} // namespace MMParametrization
} // namespace Scine