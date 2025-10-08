/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaffParameters.h"
#include "../Interactions/Electrostatic.h"
#include "../MMExceptions.h"

namespace Scine {
namespace MolecularMechanics {

std::vector<Dihedral> GaffParameters::resolveMMDihedrals(const std::string& t1, const std::string& t2,
                                                         const std::string& t3, const std::string& t4) const {
  std::vector<Dihedral> dihedralsForGivenAtomTypes;

  auto dihedralType = DihedralType(t1, t2, t3, t4);
  bool fullySpecifiedExists = (dihedrals_.find(dihedralType) != dihedrals_.end());

  // If fully specified dihedral doesn't exist, take the one with any end atom types
  if (!fullySpecifiedExists)
    dihedralType = DihedralType("X", t2, t3, "X");

  bool noDihedralTermExists = true;
  for (auto it = dihedrals_.equal_range(dihedralType).first; it != dihedrals_.equal_range(dihedralType).second; ++it) {
    noDihedralTermExists = false;
    // Add only if barrier is non-zero
    if (!(it->second.isZero()))
      dihedralsForGivenAtomTypes.push_back(it->second.toMMDihedral());
  }

  if (noDihedralTermExists)
    throw MMDihedralParametersNotAvailableException(t1, t2, t3, t4);

  return dihedralsForGivenAtomTypes;
}

std::vector<Dihedral> GaffParameters::getMMDihedrals(const std::string& t1, const std::string& t2,
                                                     const std::string& t3, const std::string& t4) const {
  std::vector<Dihedral> result;
  try {
    result = resolveMMDihedrals(t1, t2, t3, t4);
  }
  catch (const MMDihedralParametersNotAvailableException& e) {
    result = resolveMMDihedrals(getAtomClassFromType(t1), getAtomClassFromType(t2), getAtomClassFromType(t3),
                                getAtomClassFromType(t4));
  }
  return result;
}

std::vector<Dihedral> GaffParameters::getMMImproperDihedrals(const std::string& central, const std::string& t2,
                                                             const std::string& t3, const std::string& t4) const {
  std::vector<Dihedral> result = resolveMMImproperDihedrals(central, t2, t3, t4);
  if (result.empty() &&
      (type2AtomClass_.find(central) != type2AtomClass_.end() && type2AtomClass_.find(t2) != type2AtomClass_.end() &&
       type2AtomClass_.find(t3) != type2AtomClass_.end() && type2AtomClass_.find(t4) != type2AtomClass_.end())) {
    std::vector<Dihedral> classResult = resolveMMImproperDihedrals(getAtomClassFromType(central), getAtomClassFromType(t2),
                                                                   getAtomClassFromType(t3), getAtomClassFromType(t4));
    result.insert(result.end(), classResult.begin(), classResult.end());
  }
  return result;
}

std::vector<Dihedral> GaffParameters::resolveMMImproperDihedrals(const std::string& central, const std::string& t2,
                                                                 const std::string& t3, const std::string& t4) const {
  std::vector<Dihedral> improperDihedralsForGivenAtomTypes;

  auto improperDihedralType = ImproperDihedralType(central, t2, t3, t4);

  // If the specified type does not exist, look for others
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, "X", t3, t4);
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, t2, "X", t4);
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, t2, t3, "X");
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, "X", "X", t4);
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, "X", t3, "X");
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, t2, "X", "X");
  if (improperDihedrals_.find(improperDihedralType) == improperDihedrals_.end())
    improperDihedralType = ImproperDihedralType(central, "X", "X", "X");

  for (auto it = improperDihedrals_.equal_range(improperDihedralType).first;
       it != improperDihedrals_.equal_range(improperDihedralType).second; ++it) {
    improperDihedralsForGivenAtomTypes.push_back(it->second.toMMDihedral());
  }

  // NB: no error thrown if empty

  return improperDihedralsForGivenAtomTypes;
}

const LennardJonesPairParameters& GaffParameters::getMMLennardJones(const std::string& t1, const std::string& t2) const {
  const auto it1 = lennardJonesParameters_.find(t1);
  const auto it2 = lennardJonesParameters_.find(t2);
  if (it1 == lennardJonesParameters_.end() || it2 == lennardJonesParameters_.end()) {
    auto c1 = this->getAtomClassFromType(t1);
    auto c2 = this->getAtomClassFromType(t2);
    const auto it1C = lennardJonesParameters_.find(c1);
    const auto it2C = lennardJonesParameters_.find(c2);
    if (it1C == lennardJonesParameters_.end() || it2C == lennardJonesParameters_.end()) {
      throw MMLjParametersNotAvailableException(t1 + t2);
    }
    return *lennardJonesPairs_[it1C->second.getAtomTypeIndex()][it2C->second.getAtomTypeIndex()];
  }
  return *lennardJonesPairs_[it1->second.getAtomTypeIndex()][it2->second.getAtomTypeIndex()];
}

const std::string& GaffParameters::getAtomClassFromType(const std::string& type) const {
  if (type == "X") {
    return type;
  }
  const auto itC = type2AtomClass_.find(type);
  if (itC == type2AtomClass_.end()) {
    throw std::runtime_error("No atom class available for atom type " + type + ". The parameters are incomplete.");
  }
  return itC->second;
}

void GaffParameters::addDihedral(const DihedralType& dihedralType, DihedralParameters dihedralParameters) {
  dihedrals_.emplace(dihedralType, dihedralParameters);
}

void GaffParameters::addImproperDihedral(const ImproperDihedralType& improperDihedralType,
                                         DihedralParameters improperDihedralParameters) {
  improperDihedrals_.emplace(improperDihedralType, improperDihedralParameters);
}

void GaffParameters::addLennardJones(const std::string& atomType,
                                     Scine::MolecularMechanics::LennardJonesParameters lennardJonesParameters) {
  unsigned int nextIndex = lennardJonesParameters_.size();
  const bool alreadyTabulated = lennardJonesParameters_.find(atomType) != lennardJonesParameters_.end();
  lennardJonesParameters_.emplace(atomType, lennardJonesParameters);
  /*
   * Precalculate all Lennard-Jones pairs and store them in a matrix. Note that the map holding the Lennard-
   * Jones parameters may change its order if new elements are inserted. Therefore, we cache the index
   * corresponding to the Lennard-Jones parameter pairs in the matrix in the Lennard-Jones parameter object.
   */
  if (!alreadyTabulated) {
    // Add a new column and row to the parameter pair matrix.
    lennardJonesParameters_.find(atomType)->second.setAtomTypeIndex(nextIndex);
    for (auto& pairSet : lennardJonesPairs_) {
      pairSet.emplace_back(nullptr);
    }
    lennardJonesPairs_.emplace_back(nextIndex + 1, nullptr);
    // Fill the new entries.
    for (const auto& otherParams : lennardJonesParameters_) {
      unsigned int i = otherParams.second.getAtomTypeIndex();
      auto pair = std::make_shared<LennardJonesPairParameters>(otherParams.second, lennardJonesParameters);
      lennardJonesPairs_[i][nextIndex] = pair;
      lennardJonesPairs_[nextIndex][i] = pair;
    }
  }
}

bool GaffParameters::empty() const {
  return lennardJonesPairs_.empty() && dihedrals_.empty() && improperDihedrals_.empty() && bonds_.empty() && angles_.empty();
}

GaffParameters operator+(GaffParameters lhs, const GaffParameters& rhs) {
  lhs += rhs;
  return lhs;
}

GaffParameters GaffParameters::operator+=(const GaffParameters& other) {
  if (!lennardJonesPairs_.empty() && !other.getLennardJonesParameters().empty()) {
    bool inconsistentNonbonded14Scaling =
        abs(other.scalingFactorForElectrostaticOneFourTerms_ - this->scalingFactorForElectrostaticOneFourTerms_) > 1e-6 ||
        abs(other.scalingFactorForLennardJonesOneFourTerms_ - this->scalingFactorForLennardJonesOneFourTerms_) > 1e-6;
    if (inconsistentNonbonded14Scaling) {
      throw std::runtime_error("Inconsistent 1-4 non-bonded interaction scaling when combining parameter sets.");
    }
  }
  if (!other.getLennardJonesParameters().empty()) {
    this->scalingFactorForLennardJonesOneFourTerms_ = other.scalingFactorForLennardJonesOneFourTerms_;
    this->scalingFactorForElectrostaticOneFourTerms_ = other.scalingFactorForElectrostaticOneFourTerms_;
  }

  for (const auto& dihedra : other.dihedrals_) {
    this->addDihedral(dihedra.first, dihedra.second);
  }
  for (const auto& imporper : other.improperDihedrals_) {
    this->addImproperDihedral(imporper.first, imporper.second);
  }
  for (const auto& lennardJones : other.lennardJonesParameters_) {
    this->addLennardJones(lennardJones.first, lennardJones.second);
  }
  for (const auto& angle : other.angles_) {
    this->addAngle(angle.first, angle.second);
  }
  for (const auto& bond : other.bonds_) {
    this->addBond(bond.first, bond.second);
  }
  type2AtomClass_.insert(other.type2AtomClass_.begin(), other.type2AtomClass_.end());
  return *this;
}
const std::unordered_map<std::string, LennardJonesParameters>& GaffParameters::getLennardJonesParameters() const {
  return lennardJonesParameters_;
}
void GaffParameters::set14ElectrostaticScaling(double scaling) {
  scalingFactorForElectrostaticOneFourTerms_ = scaling;
}
void GaffParameters::set14LennardJonesScaling(double scaling) {
  scalingFactorForLennardJonesOneFourTerms_ = scaling;
}
double GaffParameters::get14ElectrostaticScaling() const {
  return scalingFactorForElectrostaticOneFourTerms_;
}
double GaffParameters::get14LennardJonesScaling() const {
  return scalingFactorForLennardJonesOneFourTerms_;
}
void GaffParameters::setType2AtomClass(const std::map<std::string, std::string>& map) {
  type2AtomClass_ = map;
}
Angle GaffParameters::getMMAngle(std::string t1, std::string t2, std::string t3) const {
  auto angle_ptr = angles_.find(AngleType(t1, t2, t3));
  if (angle_ptr == angles_.end()) {
    angle_ptr = angles_.find(AngleType(getAtomClassFromType(t1), getAtomClassFromType(t2), getAtomClassFromType(t3)));
  }
  if (angle_ptr == angles_.end())
    throw MMAngleParametersNotAvailableException(t1, t2, t3);
  return angle_ptr->second.toMMAngle();
}
Bond GaffParameters::getMMBond(std::string t1, std::string t2) const {
  auto bond_ptr = bonds_.find(BondType(t1, t2));
  if (bond_ptr == bonds_.end()) {
    bond_ptr = bonds_.find(BondType(getAtomClassFromType(t1), getAtomClassFromType(t2)));
  }
  if (bond_ptr == bonds_.end()) {
    throw MMBondParametersNotAvailableException(t1, t2);
  }
  return bond_ptr->second.toMMBond();
}

} // namespace MolecularMechanics
} // namespace Scine
