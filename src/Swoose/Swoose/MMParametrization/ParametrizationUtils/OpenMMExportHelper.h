/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_OPENMMEXPORTHELPER_H
#define MMPARAMETRIZATION_OPENMMEXPORTHELPER_H

#include <Swoose/MMParametrization/ParametrizationData.h>
#include <Swoose/MolecularMechanics/Interactions/HydrogenBond.h>
#include <Swoose/MolecularMechanics/SFAM/SfamPotentialTermsGenerator.h>
#include <memory>

namespace Scine {

namespace Utils {
class Settings;
} // namespace Utils

namespace MolecularMechanics {
class SfamParameters;
} // namespace MolecularMechanics

namespace MMParametrization {
struct ParametrizationData;

class OpenMMExportHelper {
 public:
  /**
   * @brief Constructor.
   */
  OpenMMExportHelper(std::shared_ptr<Utils::Settings> settings, ParametrizationData& data);
  /**
   * @brief Exports the tolopology file for OpenMM.
   */
  void exportTopology();
  /**
   * @brief Exports the SFAM parameters and functional form in the XML file format.
   */
  void exportSfamForceField();

 private:
  void writeHeader(std::ofstream& xmlFile);
  void writeAtomTypes(std::ofstream& xmlFile);
  void writeResidues(std::ofstream& xmlFile);
  void writeHarmonicBondTerms(std::ofstream& xmlFile);
  void writeAngleTerms(std::ofstream& xmlFile);
  void writeDihedralTerms(std::ofstream& xmlFile);
  void writeImproperDihedralTerms(std::ofstream& xmlFile);
  void writeCoulombTerm(std::ofstream& xmlFile);
  void writeNonbondedTerms(std::ofstream& xmlFile);
  void writeHydrogenBondTerms(std::ofstream& xmlFile);
  // calculate the cutoff Radii R_0 from Grimme's D3 model from the C6 and C8 matrix.
  Eigen::MatrixXf calculateR0Matrix();
  // check if an atom is bonded to a hydrogen atom and returns the corresponding H-index.
  std::pair<bool, int> isBondedToHydrogenAtom(int atomIndex);
  // for the definition oh hydrogen-bonded groups, check:
  // Brunken et al, JCTC, 2020, 16, 1646-1665.
  void detectDonorAndAcceptorGroups();
  // Get the unique atom types from the parametrization data and initialize further attributes.
  void initializeUniqueAtomTypes();

  // The settings
  std::shared_ptr<Utils::Settings> settings_;
  MMParametrization::ParametrizationData& data_;
  MolecularMechanics::SfamParameters parameters_;
  const unsigned int nAtoms_;
  std::vector<std::string> uniqueAtomTypes_;
  std::vector<MolecularMechanics::HydrogenBondHelper::Donor> donors_;
  std::vector<MolecularMechanics::HydrogenBondHelper::Acceptor> acceptors_;
  std::vector<double> effectiveChargesHolder_;
  std::vector<Utils::ElementType> elementTypes_;
  const std::string xmlFileName_;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_OPENMMEXPORTHELPER_H
