/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "QmmmHelpers.h"
#include "QmmmCalculatorSettings.h"
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <fstream>
#include <set>

namespace Scine {
namespace Qmmm {
namespace QmmmHelpers {

std::vector<double> writeChargesAndPositionsAsList(const Utils::AtomCollection& structure,
                                                   const ChargeRedistributionResult& chargeRedistributionResult,
                                                   const std::vector<int>& listOfQmAtoms) {
  if (!chargeRedistributionResult.atomicCharges.size() && !chargeRedistributionResult.auxiliaryCharges.size()) {
    return {};
  }
  std::vector<double> chargesAndPositions;
  const Utils::ElementTypeCollection& elements = structure.getElements();
  const Utils::PositionCollection& positions = structure.getPositions();
  for (int i = 0; i < positions.rows(); ++i) {
    if (std::find(listOfQmAtoms.begin(), listOfQmAtoms.end(), i) == listOfQmAtoms.end()) {
      const Utils::Position& pos = positions.row(i);
      const double charge = chargeRedistributionResult.atomicCharges[i];
      const auto atomicNumber = static_cast<double>(Utils::ElementInfo::Z(elements[i]));
      chargesAndPositions.push_back(charge);
      chargesAndPositions.push_back(atomicNumber);
      for (int j = 0; j < 3; ++j) {
        chargesAndPositions.push_back(pos[j]);
      }
    }
  }
  for (int i = 0; i < static_cast<int>(chargeRedistributionResult.auxiliaryCharges.size()); ++i) {
    const Utils::Position& pos = chargeRedistributionResult.positionsOfAuxiliaryCharges.row(i);
    const double charge = chargeRedistributionResult.auxiliaryCharges[i];
    chargesAndPositions.push_back(charge);
    chargesAndPositions.push_back(0.0); // 0.0 as atomic number for auxiliary charges
    for (int j = 0; j < 3; ++j) {
      chargesAndPositions.push_back(pos[j]);
    }
  }
  return chargesAndPositions;
}

void writePointChargesFile(const Utils::PositionCollection& positions,
                           const ChargeRedistributionResult& chargeRedistributionResult,
                           const std::vector<int>& listOfQmAtoms, const std::string& filename, bool writeTurbomoleFormat) {
  if (!chargeRedistributionResult.atomicCharges.size() && !chargeRedistributionResult.auxiliaryCharges.size()) {
    return;
  }
  std::ofstream pcFile(filename);
  assert(chargeRedistributionResult.atomicCharges.size() == positions.rows());
  assert(chargeRedistributionResult.auxiliaryCharges.size() == chargeRedistributionResult.positionsOfAuxiliaryCharges.rows());

  if (!writeTurbomoleFormat) {
    pcFile << chargeRedistributionResult.atomicCharges.size() - listOfQmAtoms.size() +
                  chargeRedistributionResult.auxiliaryCharges.size()
           << "\n";
  }

  // First write the atomic charges
  for (int i = 0; i < positions.rows(); ++i) {
    // Write charge if this atom is NOT a QM atom.
    double charge = chargeRedistributionResult.atomicCharges.at(i);
    // TODO: This is a quick fix to ensure that no point charges are written which are
    // so small that Turbomole doesn't consider them. A cleaner solution would be
    // desirable.
    // Note that we set the charge to 2e-6 such that all point charges are properly
    // read in: our Turbomole output parser checks for a value of 1e-6 to determine
    // the number of point charges.
    if (writeTurbomoleFormat && std::fabs(charge) < 2e-6) {
      charge = 2e-6;
    }
    if (std::find(listOfQmAtoms.begin(), listOfQmAtoms.end(), i) == listOfQmAtoms.end()) {
      if (writeTurbomoleFormat) {
        pcFile << positions.row(i) << " ";
        pcFile << charge << "\n";
      }
      else {
        pcFile << charge << " ";
        pcFile << positions.row(i) * Utils::Constants::angstrom_per_bohr << "\n";
      }
    }
  }

  // Then write the additional auxiliary charges if there are any
  for (int j = 0; j < int(chargeRedistributionResult.auxiliaryCharges.size()); ++j) {
    if (writeTurbomoleFormat) {
      pcFile << chargeRedistributionResult.positionsOfAuxiliaryCharges.row(j) * Utils::Constants::angstrom_per_bohr << " ";
      pcFile << chargeRedistributionResult.auxiliaryCharges.at(j) << "\n";
    }
    else {
      pcFile << chargeRedistributionResult.auxiliaryCharges.at(j) << " ";
      pcFile << chargeRedistributionResult.positionsOfAuxiliaryCharges.row(j) * Utils::Constants::angstrom_per_bohr << "\n";
    }
  }
}

Utils::AtomCollection createQmRegion(const std::vector<int>& listOfQmAtoms, const Utils::AtomCollection& structure,
                                     const std::vector<std::list<int>>& listsOfNeighbors,
                                     const std::string& xyzFilename, std::vector<int>& mmBoundaryAtoms) {
  Utils::AtomCollection qmRegion;
  for (const auto& qmAtomIndex : listOfQmAtoms) {
    qmRegion.push_back(structure.at(qmAtomIndex));
  }

  // Add the link atoms to preliminary QM region
  addAllLinkAtoms(qmRegion, structure, listsOfNeighbors, listOfQmAtoms, mmBoundaryAtoms);

  // Write QM region to XYZ file format if desired
  if (!xyzFilename.empty())
    Utils::ChemicalFileHandler::write(xyzFilename, qmRegion);

  return qmRegion;
}

void checkValidityOfQmRegion(const std::vector<int>& listOfQmAtoms, const Utils::AtomCollection& structure) {
  for (const auto& i : listOfQmAtoms) {
    if ((i < 0) || (i >= structure.size()))
      throw std::runtime_error(
          "The selected QM region is not valid, because at least one given atom index was negative or too large.");
  }

  const std::set<int> s(listOfQmAtoms.begin(), listOfQmAtoms.end());
  if (s.size() != listOfQmAtoms.size())
    throw std::runtime_error("The list of QM atoms contains duplicates!");

  // TODO: Add more complicated checks! (QM region should be one (or more) continuous graph(s))
}

void addAllLinkAtoms(Utils::AtomCollection& qmRegion, const Utils::AtomCollection& fullStructure,
                     const std::vector<std::list<int>>& listsOfNeighbors, const std::vector<int>& listOfQmAtoms,
                     std::vector<int>& mmBoundaryAtoms) {
  const int numQmAtomsWithoutLinkAtoms = qmRegion.size();
  // Iterate over all QM atoms
  for (int i = 0; i < numQmAtomsWithoutLinkAtoms; ++i) {
    const int qmAtomIndex = listOfQmAtoms.at(i);
    auto neighbors = listsOfNeighbors.at(qmAtomIndex);
    // Iterate over all atoms bonded to the QM atom
    for (const auto& neighbor : neighbors) {
      // Check whether this neighbor is NOT a QM atom
      if (std::find(listOfQmAtoms.begin(), listOfQmAtoms.end(), neighbor) == listOfQmAtoms.end()) {
        addOneLinkAtom(qmRegion, fullStructure.at(qmAtomIndex), fullStructure.at(neighbor));
        mmBoundaryAtoms.push_back(neighbor);
      }
    }
  }
}

void addOneLinkAtom(Utils::AtomCollection& qmRegion, const Utils::Atom& qmAtom, const Utils::Atom& mmAtom) {
  const Utils::ElementType linkAtomElement = Utils::ElementType::H;

  auto covalentRadiusQm = Utils::ElementInfo::covalentRadius(qmAtom.getElementType());
  auto covalentRadiusLink = Utils::ElementInfo::covalentRadius(linkAtomElement);

  const Eigen::RowVector3d bondVector = mmAtom.getPosition() - qmAtom.getPosition();
  const Eigen::RowVector3d scaledBondVector = ((covalentRadiusQm + covalentRadiusLink) / bondVector.norm()) * bondVector;
  const Eigen::RowVector3d linkAtomPosition = qmAtom.getPosition() + scaledBondVector;

  const Utils::Atom linkAtom(linkAtomElement, linkAtomPosition);
  qmRegion.push_back(linkAtom);
}

ChargeRedistributionResult getRedistributedCharges(std::vector<double> charges, const Utils::PositionCollection& positions,
                                                   const std::vector<int>& mmBoundaryAtoms,
                                                   const std::vector<std::list<int>>& listsOfNeighbors,
                                                   const std::vector<int>& listOfQmAtoms, const std::string& scheme) {
  ChargeRedistributionResult result;
  // Iterate over all MM boundary atoms
  for (const auto& atomIndex : mmBoundaryAtoms) {
    std::vector<int> nonQmNeighbors;
    // Iterate over all covalently bonded neighbors and add them to the vector above if they are not QM atoms
    for (const auto& neighbor : listsOfNeighbors.at(atomIndex)) {
      if (std::find(listOfQmAtoms.begin(), listOfQmAtoms.end(), neighbor) == listOfQmAtoms.end())
        nonQmNeighbors.push_back(neighbor);
    }
    // Redistribute the charge
    if (scheme == SwooseUtilities::OptionNames::redistributedChargeOption) {
      // Just add redistribute all of the charge of the MM boundary atom to its non-QM neighbors equally:
      for (const auto& nonQmNeighbor : nonQmNeighbors) {
        charges.at(nonQmNeighbor) += charges.at(atomIndex) / nonQmNeighbors.size();
      }
      // Set MM boundary charge to zero
      charges.at(atomIndex) = 0.0;
    }
    else if (scheme == SwooseUtilities::OptionNames::redistributedChargeAndDipolesOption) {
      // RCD scheme, according to J. Phys. Chem. A 2005, 109, 3991-4004.
      for (const auto& nonQmNeighbor : nonQmNeighbors) {
        // Subtract q/n from the neighbor
        charges.at(nonQmNeighbor) -= charges.at(atomIndex) / nonQmNeighbors.size();
        // Add 2q/n to the center of the corresponding bond vector
        result.auxiliaryCharges.push_back(2 * charges.at(atomIndex) / nonQmNeighbors.size());
        const Eigen::RowVector3d pos = 0.5 * (positions.row(atomIndex) + positions.row(nonQmNeighbor));
        result.positionsOfAuxiliaryCharges.conservativeResize(result.positionsOfAuxiliaryCharges.rows() + 1, 3);
        result.positionsOfAuxiliaryCharges.row(result.positionsOfAuxiliaryCharges.rows() - 1) = pos;
      }
      // Set MM boundary charge to zero
      charges.at(atomIndex) = 0.0;
    }
    else {
      throw std::runtime_error("This charge redistribution option is not implemented!");
    }
  }

  result.atomicCharges = std::move(charges);
  return result;
}

} // namespace QmmmHelpers
} // namespace Qmmm
} // namespace Scine
