/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaffPotentialTermsGenerator.h"
#include "../AtomTypesHolder.h"
#include "../MMExceptions.h"
#include "../PotentialTermsHelper.h"
#include "../Topology/IndexedStructuralTopology.h"
#include "GaffParameters.h"
#include <Core/Log.h>

namespace Scine {
namespace MolecularMechanics {

GaffPotentialTermsGenerator::GaffPotentialTermsGenerator(int nAtoms, const AtomTypesHolder& atomTypes,
                                                         const IndexedStructuralTopology& topology,
                                                         const GaffParameters& parameters,
                                                         const Utils::PositionCollection& positions,
                                                         const double& nonCovalentCutoffRadius, Core::Log& log)
  : nAtoms_(nAtoms),
    atomTypesHolder_(atomTypes),
    topology_(topology),
    parameters_(parameters),
    positions_(positions),
    cutoff_(std::make_shared<double>(nonCovalentCutoffRadius)),
    log_(log) {
}

std::vector<BondedTerm> GaffPotentialTermsGenerator::getBondedTerms() {
  return PotentialTermsHelper::getBondedTerms(topology_, parameters_, atomTypesHolder_);
}

std::vector<AngleTerm> GaffPotentialTermsGenerator::getAngleTerms() {
  return PotentialTermsHelper::getAngleTerms(topology_, parameters_, atomTypesHolder_);
}

std::vector<DihedralTerm> GaffPotentialTermsGenerator::getDihedralTerms() {
  std::vector<DihedralTerm> dihedralList;

  for (const auto& dihedral : topology_.getDihedralContainer()) {
    auto atomType1 = atomTypesHolder_.getAtomType(dihedral.atom1);
    auto atomType2 = atomTypesHolder_.getAtomType(dihedral.atom2);
    auto atomType3 = atomTypesHolder_.getAtomType(dihedral.atom3);
    auto atomType4 = atomTypesHolder_.getAtomType(dihedral.atom4);
    std::vector<Dihedral> m;
    DihedralType dihedralType(atomType1, atomType2, atomType3, atomType4);
    try {
      // returns vector of dihedral parameters for given atom types, since it is
      // possible there are more than one in case of Fourier series
      m = parameters_.getMMDihedrals(atomType1, atomType2, atomType3, atomType4);
      for (auto& d : m) {
        d.setCosinePreFactor(1.0);
        DihedralTerm term(dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4, d, dihedralType);
        dihedralList.push_back(term);
      }
    }
    catch (const MMDihedralParametersNotAvailableException& e) {
      Dihedral d;
      d.setCosinePreFactor(1.0);
      DihedralTerm term(dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4, d, dihedralType);
      dihedralList.push_back(term);
    }
  }

  return dihedralList;
}

std::vector<DihedralTerm> GaffPotentialTermsGenerator::getImproperDihedralTerms() {
  std::vector<DihedralTerm> improperDihedralList;

  for (const auto& improperDihedral : topology_.getImproperDihedralContainer()) {
    auto atomTypeC = atomTypesHolder_.getAtomType(improperDihedral.centralAtom);
    auto atomType2 = atomTypesHolder_.getAtomType(improperDihedral.atom2);
    auto atomType3 = atomTypesHolder_.getAtomType(improperDihedral.atom3);
    auto atomType4 = atomTypesHolder_.getAtomType(improperDihedral.atom4);
    DihedralType improperDihedralType(atomType2, atomType3, atomTypeC, atomType4);

    std::vector<Dihedral> m;
    try {
      m = parameters_.getMMImproperDihedrals(atomTypeC, atomType2, atomType3, atomType4);
      for (auto& d : m) {
        d.setCosinePreFactor(1.0);
        DihedralTerm term(improperDihedral.atom2, improperDihedral.atom3, improperDihedral.centralAtom,
                          improperDihedral.atom4, d, improperDihedralType);
        improperDihedralList.push_back(term);
      }
    }
    catch (const MMImproperDihedralParametersNotAvailableException& e) {
      // Currently, no error is thrown by the getMMImproperDihedral() function.
    }
  }

  return improperDihedralList;
}

std::vector<LennardJonesTerm> GaffPotentialTermsGenerator::getLennardJonesTerms(bool applyCutoff) {
  std::vector<LennardJonesTerm> ljList;

  // 0 = excluded, 1 = included, -1 = scaled
  Eigen::MatrixXi exclusionType = PotentialTermsHelper::getExclusionTypeMatrix(topology_, nAtoms_);

  for (int a = 0; a < nAtoms_; ++a) {
    for (int b = 0; b < a; ++b) {
      if (applyCutoff) {
        const auto& R = (positions_.row(b) - positions_.row(a)).norm();
        if (R > *cutoff_)
          continue;
      }

      auto type = exclusionType(a, b);
      if (type == 0)
        continue;

      auto atomType1 = atomTypesHolder_.getAtomType(a);
      auto atomType2 = atomTypesHolder_.getAtomType(b);

      double factor = 1.0;
      if (type == -1)
        factor = scalingFactorForLennardJonesOneFourTerms_;
      LennardJones lj = parameters_.getMMLennardJones(atomType1, atomType2, factor);
      LennardJonesTerm term(a, b, lj, cutoff_);
      ljList.push_back(term);
    }
  }

  return ljList;
}

std::vector<ElectrostaticTerm> GaffPotentialTermsGenerator::getElectrostaticTerms(bool applyCutoff) {
  // 0 = excluded, 1 = included, -1 = scaled
  Eigen::MatrixXi exclusionType = PotentialTermsHelper::getExclusionTypeMatrix(topology_, nAtoms_);
  return PotentialTermsHelper::getElectrostaticTerms(applyCutoff, cutoff_, scalingFactorForElectrostaticOneFourTerms_,
                                                     exclusionType, positions_);
}

} // namespace MolecularMechanics
} // namespace Scine
