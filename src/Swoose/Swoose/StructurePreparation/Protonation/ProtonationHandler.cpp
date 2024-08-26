/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ProtonationHandler.h"
#include "../SpecialCaseHandler.h"
#include "../StructurePreparationHelper.h"
#include "../StructurePreparationIO.h"
#include "../StructurePreparationSettings.h"
#include "ProtonationHelper.h"
#include <Molassembler/Graph.h>
#include <Molassembler/Interpret.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondDetector.h>
#include <Utils/Geometry/Atom.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Geometry/StructuralCompletion.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/ChemicalFileFormats/OpenBabelStreamHandler.h>
#include <Utils/IO/NativeFilenames.h>
#include <Eigen/Core>
#include <boost/filesystem.hpp>
#include <boost/process.hpp>

namespace bp = boost::process;

namespace Scine {
namespace StructurePreparation {
using namespace SpecialCaseHandler;

ProtonationHandler::ProtonationHandler(StructurePreparationData& data, StructurePreparationFiles& files,
                                       std::shared_ptr<Utils::Settings> settings)
  : files_(files), data_(data), settings_(settings) {
}

ProtonationHandler::ProtonationHandler() {
}

void ProtonationHandler::setup(const Utils::AtomCollection& structure, bool performGraphAnalysis) {
  nAtoms_ = structure.size();
  nNeighbors_.resize(nAtoms_);
  Utils::BondOrderCollection bondOrders;
  bondOrders.resize(nAtoms_);
  data_.fullStructure = structure;
  bondOrders = Utils::BondDetector::detectBonds(structure);
  listsOfNeighbors_ = SwooseUtilities::TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(
      nAtoms_, bondOrders, minimalBondOrderToConsider_);

  auto result = Scine::Molassembler::Interpret::graphs(structure, bondOrders);
  data_.atomIndexMapping.resize(nAtoms_);
  for (int i = 0; i < nAtoms_; ++i) {
    data_.atomIndexMapping[result.componentMap.apply(i).component].push_back(i);
  }
  SwooseUtilities::TopologyUtils::calculateNumberOfNeighbors(nAtoms_, listsOfNeighbors_, nNeighbors_);

  if (performGraphAnalysis) {
    StructurePreparationHelper::performGraphAnalysisOnStructure(data_, structure);
  }

  protonatedProtein_.clear();
  protonationTypes_.tetrahedral.clear();
  protonationTypes_.trigonalPlanar.clear();
  protonationTypes_.linear.clear();
  protonationTypes_.pseudoTetrahedral.clear();
  hydrogenAtoms_.clear();
}

void ProtonationHandler::protonateAllAminoAcids() {
  setup(protein_);
  // Sort the hydrogens
  sortAtomTypesByProtonationCriteria();
  // protonate all sp3 hybrid. groups
  protonateTetrahedralGroups(protein_);
  // protonate all sp2 hybrid. groups
  protonateTrigonalPlanarGroups(protein_);
  // protonate all terminal O-X and S-X
  protonateTerminalSAndO(protein_);
  // protonate peptid bonds
  protonatePeptidBonds(protein_);
  // merge hydrogens and protein
  protonatedProtein_ += protein_;
  protonatedProtein_ += hydrogenAtoms_;
  Utils::ChemicalFileHandler::write(files_.protonatedProteinFile, protonatedProtein_);
}

void ProtonationHandler::protonateNonRegContainerWithExternalSoftware() {
  if (Utils::OpenBabelStreamHandler::checkForBinary()) {
    checkObabelVersion();
    if (boost::filesystem::exists(files_.nonRegContainerFile)) {
      // Convert xyz to Pdb internally because the obabel protonation is then more reliable
      std::string nonRegContainerPdb =
          Utils::NativeFilenames::combinePathSegments(files_.workingDirectory, "nonRegContainer.pdb");
      StructurePreparationIO::xyzToPdb(files_.nonRegContainerFile, nonRegContainerPdb);
      // Protonate NonRegContainer with obabel
      std::string obabelExecutableString =
          "obabel -ipdb " + nonRegContainerPdb + " -O " + files_.protonatedNonRegContainerFile + " -h -p " +
          std::to_string(settings_->getDouble(SwooseUtilities::SettingsNames::phValueOfSystem));
      bp::ipstream err;
      bp::system(obabelExecutableString, bp::std_out > files_.protonatedNonRegContainerFile, bp::std_err > err);
      if (!ProtonationHelper::openBabelSuccess(err)) {
        throw std::runtime_error(
            "Error encountered whilst protonating the nonRegContainer with OpenBabel! Please check if "
            "you installed OpenBabel correctly!");
      }
      boost::filesystem::remove(nonRegContainerPdb);
      // Now write empty atomic info file and let the user modify it.
      std::ofstream atomicInfoNonRegContainer(files_.nonRegContainerInfoFile);
    }
  }
  else {
    throw std::runtime_error("OpenBabel is not available. Please download obabel and place it in your path. ");
  }
}

void ProtonationHandler::sortAtomTypesByProtonationCriteria() {
  for (const auto& proteinAtom : data_.protein) {
    // Atom types that are always tetrahedral.
    int numNeighbors = nNeighbors_[proteinAtom.index];
    if (ProtonationHelper::isType(listOfsp3AtomTypes_, proteinAtom.atomType) && (numNeighbors < 4)) {
      protonationTypes_.tetrahedral.push_back(proteinAtom.index);
    }
    // Atom types that are always trigonal planar.
    else if (ProtonationHelper::isType(listOfsp2AtomTypes_, proteinAtom.atomType) && (numNeighbors < 3)) {
      protonationTypes_.trigonalPlanar.push_back(proteinAtom.index);
    }
    // Atom types that are always pseudo-tetrahedral.
    else if (ProtonationHelper::isType(listOfPseudoTetrahedrals_, proteinAtom.atomType) && (numNeighbors < 3)) {
      protonationTypes_.pseudoTetrahedral.push_back(proteinAtom.index);
    }
    else if (proteinAtom.atomType == "CG") {
      if (ProtonationHelper::isType(AAwithsp3CG_, proteinAtom.residueName)) {
        protonationTypes_.tetrahedral.push_back(proteinAtom.index);
      }
      else if (numNeighbors < 3) {
        protonationTypes_.trigonalPlanar.push_back(proteinAtom.index);
      }
    }
    // CD
    else if (proteinAtom.atomType == "CD") {
      if (ProtonationHelper::isType(AAwithsp3CD_, proteinAtom.residueName)) {
        protonationTypes_.tetrahedral.push_back(proteinAtom.index);
      }
      else if (numNeighbors < 3) {
        protonationTypes_.trigonalPlanar.push_back(proteinAtom.index);
      }
    }
    // CD1 and CD2
    else if (proteinAtom.atomType == "CD1" || proteinAtom.atomType == "CD2") {
      if (proteinAtom.residueName == "LEU" || proteinAtom.residueName == "ILE") {
        protonationTypes_.tetrahedral.push_back(proteinAtom.index);
      }
      else {
        if (numNeighbors < 3) {
          protonationTypes_.trigonalPlanar.push_back(proteinAtom.index);
        }
      }
    }
    // NE2 in GLN
    else if (proteinAtom.atomType == "NE2" && (proteinAtom.residueName == "GLN")) {
      protonationTypes_.trigonalPlanar.push_back(proteinAtom.index);
    }

    else if (ProtonationHelper::isType(terminalOhSh_, proteinAtom.atomType)) {
      protonationTypes_.linear.push_back(proteinAtom.index);
    }
    else if (proteinAtom.atomType == "SG" && numNeighbors == 1) {
      protonationTypes_.linear.push_back(proteinAtom.index);
    }
    else if (proteinAtom.atomType == "N") {
      PeptidBond bond;
      bond.N = proteinAtom.index;
      bool isValidPeptid = false;
      getPeptidBonds(bond.N, bond.C, bond.CA, bond.O, isValidPeptid);
      if (isValidPeptid) {
        if (proteinAtom.residueName != "PRO") {
          peptidBonds_.push_back(bond);
        }
      }
      else {
        if (settings_->getBool(SwooseUtilities::SettingsNames::chargedCandNTermini)) {
          protonationTypes_.tetrahedral.push_back(proteinAtom.index);
        }
        else {
          protonationTypes_.tetrahedral.push_back(proteinAtom.index);
        }
      }
    }
    // saturate incomplete C termini as aldehydes
    else if (proteinAtom.atomType == "C" && numNeighbors == 2) {
      protonationTypes_.trigonalPlanar.push_back(proteinAtom.index);
    }
    else if ((proteinAtom.atomType == "OD2" && proteinAtom.residueName == "ASP") ||
             (proteinAtom.atomType == "OE2" && proteinAtom.residueName == "GLU")) {
      protonationTypes_.linear.push_back(proteinAtom.index);
    }
    detectCTerminus(proteinAtom);
  }
}

void ProtonationHandler::detectCTerminus(const ProteinAtom& proteinAtom) {
  if (proteinAtom.atomType == "C") {
    int numNeighbors = nNeighbors_[proteinAtom.index];
    if (numNeighbors == 3) {
      int protonatedOxygens = 0;
      auto neighborsOfC = listsOfNeighbors_[proteinAtom.index];
      for (auto& neighbor : neighborsOfC) {
        if (data_.protein[neighbor].atomType == "O" && protonatedOxygens > 1)
          if (!settings_->getBool(SwooseUtilities::SettingsNames::chargedCandNTermini)) {
            protonationTypes_.linear.push_back(data_.protein[neighbor].index);
          }
        protonatedOxygens++;
      }
    }
    else if (numNeighbors == 2) {
      protonationTypes_.trigonalPlanar.push_back(proteinAtom.index);
    }
  }
}

void ProtonationHandler::getPeptidBonds(int indexN, int& indexC, int& indexCA, int& indexO, bool& isValidPeptid) {
  auto neighborsOfN = listsOfNeighbors_[indexN];

  isValidPeptid = false;
  bool nitrogenHasCaNeighbor = false;
  bool nitrogenHasCNeighbor = false;
  for (auto neighbor : neighborsOfN) {
    if (data_.protein[neighbor].atomType == "CA") {
      nitrogenHasCaNeighbor = true;
      indexCA = neighbor;
    }
    else if (data_.protein[neighbor].atomType == "C") {
      nitrogenHasCNeighbor = true;
      indexC = neighbor;
    }
  }

  if (nitrogenHasCaNeighbor && nitrogenHasCNeighbor) {
    isValidPeptid = true;
    auto neighborsOfC = listsOfNeighbors_[indexC];
    for (const auto neighborOfC : neighborsOfC) {
      if (protein_.getElement(neighborOfC) == Utils::ElementType::O) {
        indexO = neighborOfC;
      }
    }
  }
}

void ProtonationHandler::reprotonateCriticalAtom(Utils::AtomCollection& structure, int indexOfCriticalAtom,
                                                 const std::string& residueName) {
  bool performGraphAnalysis = false;
  setup(structure, performGraphAnalysis);
  if (residueName == "ARG" || residueName == "LYS") {
    protonationTypes_.tetrahedral.push_back(indexOfCriticalAtom);
  }
  else if (residueName == "HIS") {
    protonationTypes_.trigonalPlanar.push_back(indexOfCriticalAtom);
  }

  protonateTetrahedralGroups(structure);
  protonateTrigonalPlanarGroups(structure);
  structure += hydrogenAtoms_;
}

Utils::Atom ProtonationHandler::generatePeptidHydrogen(const Utils::AtomCollection& structure, int indexN, int indexC,
                                                       int indexO) {
  Utils::ElementType linkAtomElement = Utils::ElementType::H;

  Eigen::RowVector3d cVector = structure.getPosition(indexC);
  Eigen::RowVector3d nVector = structure.getPosition(indexN);
  Eigen::RowVector3d oVector = structure.getPosition(indexO);

  Eigen::RowVector3d coBondVector = oVector - cVector;
  Eigen::RowVector3d inverseScaledBondVector = -coBondVector / coBondVector.norm() * 2.0;
  Eigen::RowVector3d linkAtomPosition = nVector + inverseScaledBondVector;

  Utils::Atom hydrogenAtom(linkAtomElement, linkAtomPosition);
  return hydrogenAtom;
}

void ProtonationHandler::protonatePeptidBonds(const Utils::AtomCollection& structure) {
  for (auto peptid : peptidBonds_) {
    auto hydrogen = generatePeptidHydrogen(structure, peptid.N, peptid.C, peptid.O);
    hydrogenAtoms_.push_back(hydrogen);
  }
}

void ProtonationHandler::protonateTetrahedralGroups(const Utils::AtomCollection& structure) {
  for (int index = 0; index < structure.size(); ++index) {
    bool isTetrahedralAtom = ProtonationHelper::isAtomOf(protonationTypes_.tetrahedral, index);
    bool isPseudoTetrahedralAtom = ProtonationHelper::isAtomOf(protonationTypes_.pseudoTetrahedral, index);
    if (isTetrahedralAtom || isPseudoTetrahedralAtom) {
      auto pos = structure.at(index).getPosition();

      if (nNeighbors_[index] == 1) {
        auto firstBondedAtom = listsOfNeighbors_[index].begin();
        std::vector<Eigen::Vector3d> hPos(3); // v2, v3, v4;
        const Eigen::Vector3d v1 = (structure.at(*firstBondedAtom).getPosition() - pos).normalized();

        Utils::StructuralCompletion::generate3TetrahedronCornersFrom1Other(v1, hPos.at(0), hPos.at(1), hPos.at(2));
        size_t counter = 0;
        for (auto i : hPos) {
          if (counter == hPos.size() - 1 && isPseudoTetrahedralAtom) {
            break;
          }
          auto bonded_element = structure.at(index).getElementType();
          const double sumOfVdWRadii = Utils::ElementInfo::covalentRadius(bonded_element) +
                                       Utils::ElementInfo::covalentRadius(Utils::ElementType::H);
          i *= sumOfVdWRadii;
          i += pos;
          const Utils::Atom hydrogen(Utils::ElementType::H, i);
          hydrogenAtoms_.push_back(hydrogen);
          counter++;
        }
      }
      else if (nNeighbors_[index] == 2) {
        auto firstBondedAtom = listsOfNeighbors_[index].begin();
        auto secondBondedAtom = std::next(listsOfNeighbors_[index].begin());

        std::vector<Eigen::Vector3d> hPos(2);
        const Eigen::Vector3d v1 = (structure.at(*firstBondedAtom).getPosition() - pos).normalized();
        const Eigen::Vector3d v2 = (structure.at(*secondBondedAtom).getPosition() - pos).normalized();

        Utils::StructuralCompletion::generate2TetrahedronCornersFrom2Others(v1, v2, hPos.at(0), hPos.at(1));
        // TODO: Is there a reason why this counter is never incremented?
        size_t counter = 0;
        for (auto i : hPos) {
          if (counter == hPos.size() - 1 && isPseudoTetrahedralAtom) {
            break;
          }
          auto bonded_element = structure.at(index).getElementType();
          const double sumOfVdWRadii = Utils::ElementInfo::covalentRadius(bonded_element) +
                                       Utils::ElementInfo::covalentRadius(Utils::ElementType::H);
          i *= sumOfVdWRadii;
          i += pos;
          const Utils::Atom hydrogen(Utils::ElementType::H, i);

          hydrogenAtoms_.push_back(hydrogen);
        }
      }

      else if (nNeighbors_[index] == 3) {
        auto firstBondedAtom = listsOfNeighbors_[index].begin();
        auto secondBondedAtom = firstBondedAtom++;
        auto thirdBondedAtom = firstBondedAtom++;
        Eigen::Vector3d v4;
        Eigen::Vector3d v1 = (structure.at(*firstBondedAtom).getPosition() - pos).normalized();
        Eigen::Vector3d v2 = (structure.at(*secondBondedAtom).getPosition() - pos).normalized();
        Eigen::Vector3d v3 = (structure.at(*thirdBondedAtom).getPosition() - pos).normalized();

        Utils::StructuralCompletion::generate1TetrahedronCornerFrom3Others(v1, v2, v3, v4);
        auto bonded_element = structure.at(index).getElementType();
        double sumOfVdWRadii = Utils::ElementInfo::covalentRadius(bonded_element) +
                               Utils::ElementInfo::covalentRadius(Utils::ElementType::H);
        v4 *= sumOfVdWRadii;
        v4 += pos;
        Utils::Atom hydrogen(Utils::ElementType::H, v4);
        hydrogenAtoms_.push_back(hydrogen);
      }
      else {
        throw std::runtime_error("Protonation failed for tetrahedral group. " + std::to_string(nNeighbors_[index]));
      }
    }
  }
}

void ProtonationHandler::protonateTrigonalPlanarGroups(const Utils::AtomCollection& structure) {
  for (int index = 0; index < structure.size(); ++index) {
    if (ProtonationHelper::isAtomOf(protonationTypes_.trigonalPlanar, index)) {
      auto pos = structure.at(index).getPosition();

      if (nNeighbors_[index] == 1) {
        auto firstBondedAtom = listsOfNeighbors_[index].begin();
        std::vector<Eigen::Vector3d> hPos(2);
        Eigen::Vector3d v1 = (structure.at(*firstBondedAtom).getPosition() - pos).normalized();
        Utils::StructuralCompletion::generate2TriangleCornersFrom1Other(v1, hPos.at(0), hPos.at(1));
        for (auto i : hPos) {
          i.normalize();
          auto bonded_element = structure.at(index).getElementType();
          double sumOfVdWRadii = Utils::ElementInfo::covalentRadius(bonded_element) +
                                 Utils::ElementInfo::covalentRadius(Utils::ElementType::H);
          i *= sumOfVdWRadii;
          i += pos;
          Utils::Atom hydrogen(Utils::ElementType::H, i);
          hydrogenAtoms_.push_back(hydrogen);
        }
      }

      else if (nNeighbors_[index] == 2) {
        auto firstBondedAtom = listsOfNeighbors_[index].begin();
        auto secondBondedAtom = firstBondedAtom++;

        Eigen::Vector3d v3;
        Eigen::Vector3d v1 = (structure.at(*firstBondedAtom).getPosition() - pos).normalized();
        Eigen::Vector3d v2 = (structure.at(*secondBondedAtom).getPosition() - pos).normalized();

        Utils::StructuralCompletion::generate1TriangleCornerFrom2Others(v1, v2, v3);
        v3.normalize();
        auto bonded_element = structure.at(index).getElementType();
        double sumOfVdWRadii = Utils::ElementInfo::covalentRadius(bonded_element) +
                               Utils::ElementInfo::covalentRadius(Utils::ElementType::H);
        v3 *= sumOfVdWRadii;
        v3 += pos;
        Utils::Atom hydrogen(Utils::ElementType::H, v3);

        hydrogenAtoms_.push_back(hydrogen);
      }
      else {
        throw std::runtime_error("Protonation failed for a trigonal planar group. ");
      }
    }
  }
}

void ProtonationHandler::protonateTerminalSAndO(const Utils::AtomCollection& structure) {
  for (int index = 0; index < structure.size(); ++index) {
    if (ProtonationHelper::isAtomOf(protonationTypes_.linear, index)) {
      auto pos = structure.at(index).getPosition();

      if (nNeighbors_[index] == 1) {
        auto firstBondedAtom = listsOfNeighbors_[index].begin();
        std::vector<Eigen::Vector3d> hPos(3);
        Eigen::Vector3d v1 = (structure.at(*firstBondedAtom).getPosition() - pos).normalized();
        Utils::StructuralCompletion::generate3TetrahedronCornersFrom1Other(v1, hPos.at(0), hPos.at(1), hPos.at(2));
        auto bonded_element = structure.at(index).getElementType();
        double sumOfVdWRadii = Utils::ElementInfo::covalentRadius(bonded_element) +
                               Utils::ElementInfo::covalentRadius(Utils::ElementType::H);
        hPos.at(0) *= sumOfVdWRadii;
        hPos.at(0) += pos;
        Utils::Atom hydrogen(Utils::ElementType::H, hPos.at(0));
        hydrogenAtoms_.push_back(hydrogen);
      }
    }
  }
}

void ProtonationHandler::setProteinStructure(Utils::AtomCollection& structure) {
  protein_ = std::move(structure);
}

void ProtonationHandler::checkObabelVersion() {
  std::string callString = "obabel -V";
  boost::process::pstream readVersionStream;
  boost::process::child readVersionChild(callString, boost::process::std_out > readVersionStream);
  readVersionChild.wait();

  std::string versionString;
  std::getline(readVersionStream, versionString);

  // Note that, for some reason, also version 3.1.1 has a version string of "3.1.0"
  if (versionString.find("3.0.0") == std::string::npos && versionString.find("3.1.0") == std::string::npos &&
      versionString.find("3.1.1") == std::string::npos) {
    throw std::runtime_error("Only versions 3.0.0, 3.1.0, and 3.1.1 of OpenBabel are supported.");
  }
}

} // namespace StructurePreparation
} // namespace Scine
