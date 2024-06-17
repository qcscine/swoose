/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SfamParameters.h"
#include "../AtomTypesHolder.h"
#include "../Interactions/Electrostatic.h"
#include "../MMExceptions.h"
#include <Utils/Geometry/AtomCollection.h>

namespace Scine {
namespace MolecularMechanics {

bool SfamParameters::sanityCheck(const AtomTypesHolder& atomTypes) const {
  int numDistinctAtomTypes = evaluateNumDistinctAtomTypes(atomTypes);
  if (c6Matrix_.rows() != numDistinctAtomTypes) {
    return false;
  }
  if (c6Matrix_.cols() != numDistinctAtomTypes) {
    return false;
  }
  if (int(c6IndicesMap_.size()) != numDistinctAtomTypes) {
    return false;
  }
  if (nonCovalentParameters_.size() != 5) {
    return false;
  }
  return true;
}

std::vector<double> SfamParameters::getChargesForEachAtom(const AtomTypesHolder& atomTypes) const {
  std::vector<double> charges;
  charges.reserve(atomTypes.size());
  for (int i = 0; i < atomTypes.size(); ++i) {
    try {
      double charge = charges_.at(atomTypes.getAtomType(i));
      charges.push_back(charge);
    }
    catch (const std::out_of_range& e) {
      throw std::runtime_error("There is no atomic charge parameter available for atom " + std::to_string(i) + " " +
                               atomTypes.getAtomType(i));
    }
  }
  return charges;
}

// c6 coefficients are stored as floats for efficiency
void SfamParameters::prepareC6Matrix(const AtomTypesHolder& atomTypes) {
  c6IndicesMap_.clear();
  // Add atom types from vector to map. Afterwards, they'll be ordered.
  for (int i = 0; i < atomTypes.size(); ++i)
    c6IndicesMap_[atomTypes.getAtomType(i)] = 0;
  // Adjust the indices
  int idx = 0;
  for (auto& element : c6IndicesMap_) {
    element.second = idx;
    idx++;
  }
  c6Matrix_.resize(idx, idx);
  c6Matrix_.setZero();
}

void SfamParameters::resetC6Matrix() {
  c6Matrix_.resize(0, 0);
  c6IndicesMap_.clear();
}

// c6 coefficients are stored as floats for efficiency
float SfamParameters::getC6(const std::string& a, const std::string& b) const {
  if (c6IndicesMap_.find(a) == c6IndicesMap_.end() || c6IndicesMap_.find(b) == c6IndicesMap_.end())
    throw std::runtime_error("C6 coefficient was requested for pair " + a + " - " + b +
                             ". At least one of the two atom types are unknown.");
  return c6Matrix_(c6IndicesMap_.at(a), c6IndicesMap_.at(b));
}

float SfamParameters::getC6(int indexOfAtomTypeA, int indexOfAtomTypeB) const {
  return c6Matrix_(indexOfAtomTypeA, indexOfAtomTypeB);
}

void SfamParameters::setC6(const std::string& a, const std::string& b, float c6) {
  if (c6IndicesMap_.find(a) == c6IndicesMap_.end() || c6IndicesMap_.find(b) == c6IndicesMap_.end())
    throw std::runtime_error("C6 coefficient was set for pair " + a + " - " + b +
                             ". At least one of the two atom types are unknown.");
  c6Matrix_(c6IndicesMap_[a], c6IndicesMap_[b]) = c6;
  c6Matrix_(c6IndicesMap_[b], c6IndicesMap_[a]) = c6;
}

void SfamParameters::setC6(int indexOfAtomTypeA, int indexOfAtomTypeB, float c6) {
  c6Matrix_(indexOfAtomTypeA, indexOfAtomTypeB) = c6;
  c6Matrix_(indexOfAtomTypeB, indexOfAtomTypeA) = c6;
}

std::vector<double> SfamParameters::getNonCovalentParameters() const {
  return nonCovalentParameters_;
}

void SfamParameters::setNonCovalentParameters(std::vector<double> nonCovalentParameters) {
  nonCovalentParameters_ = std::move(nonCovalentParameters);
}

std::vector<Dihedral> SfamParameters::getMMDihedrals(std::string t1, std::string t2, std::string t3, std::string t4) const {
  std::vector<Dihedral> dihedralsForGivenAtoms;

  auto dihedralType = DihedralType("X", t2, t3, "X");

  bool noDihedralTermExists = true;
  for (auto it = dihedrals_.equal_range(dihedralType).first; it != dihedrals_.equal_range(dihedralType).second; ++it) {
    noDihedralTermExists = false;
    // Add only if barrier is non-zero
    if (!(it->second.isZero())) {
      dihedralsForGivenAtoms.push_back(it->second.toMMDihedral());
    }
  }

  if (noDihedralTermExists)
    throw MMDihedralParametersNotAvailableException(t1, t2, t3, t4);

  return dihedralsForGivenAtoms;
}

ImproperDihedral SfamParameters::getMMImproperDihedral(std::string central, std::string t2, std::string t3, std::string t4) const {
  auto improper_dihedral_ptr = improperDihedrals_.find(ImproperDihedralType(central, t2, t3, t4));
  if (improper_dihedral_ptr == improperDihedrals_.end())
    throw MMImproperDihedralParametersNotAvailableException(central, t2, t3, t4);
  else
    return improper_dihedral_ptr->second.toMMImproperDihedral();
}

void SfamParameters::addCharge(std::string atomType, double charge) {
  charges_.emplace(atomType, charge);
}

void SfamParameters::addDihedral(DihedralType dihedralType, DihedralParameters dihedralParameters) {
  dihedrals_.emplace(dihedralType, dihedralParameters);
}

void SfamParameters::addImproperDihedral(ImproperDihedralType improperDihedralType,
                                         ImproperDihedralParameters improperDihedralParameters) {
  improperDihedrals_.emplace(improperDihedralType, improperDihedralParameters);
}

std::map<BondType, BondParameters>& SfamParameters::getBonds() {
  return bonds_;
}

const std::map<BondType, BondParameters>& SfamParameters::getBonds() const {
  return bonds_;
}

std::map<AngleType, AngleParameters>& SfamParameters::getAngles() {
  return angles_;
}

const std::map<AngleType, AngleParameters>& SfamParameters::getAngles() const {
  return angles_;
}

std::map<DihedralType, DihedralParameters>& SfamParameters::getDihedrals() {
  return dihedrals_;
}

const std::map<DihedralType, DihedralParameters>& SfamParameters::getDihedrals() const {
  return dihedrals_;
}

std::map<ImproperDihedralType, ImproperDihedralParameters>& SfamParameters::getImproperDihedrals() {
  return improperDihedrals_;
}

const std::map<ImproperDihedralType, ImproperDihedralParameters>& SfamParameters::getImproperDihedrals() const {
  return improperDihedrals_;
}

std::map<std::string, double>& SfamParameters::getCharges() {
  return charges_;
}

const std::map<std::string, double>& SfamParameters::getCharges() const {
  return charges_;
}

const std::map<std::string, int>& SfamParameters::getC6IndicesMap() const {
  return c6IndicesMap_;
}

int SfamParameters::evaluateNumDistinctAtomTypes(const AtomTypesHolder& atomTypes) const {
  std::vector<std::string> distinctAtomTypes;
  for (int i = 0; i < atomTypes.size(); ++i) {
    if (std::find(distinctAtomTypes.begin(), distinctAtomTypes.end(), atomTypes.getAtomType(i)) == distinctAtomTypes.end()) {
      distinctAtomTypes.push_back(atomTypes.getAtomType(i));
    }
  }
  return distinctAtomTypes.size();
}

} // namespace MolecularMechanics
} // namespace Scine
