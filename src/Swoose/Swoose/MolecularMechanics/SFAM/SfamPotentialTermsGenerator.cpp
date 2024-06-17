/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SfamPotentialTermsGenerator.h"
#include "../AtomTypesHolder.h"
#include "../MMExceptions.h"
#include "../PotentialTermsHelper.h"
#include "../Topology/IndexedStructuralTopology.h"
#include "SfamParameters.h"
#include <Core/Log.h>

namespace Scine {
namespace MolecularMechanics {

SfamPotentialTermsGenerator::SfamPotentialTermsGenerator(int nAtoms, const AtomTypesHolder& atomTypes,
                                                         const IndexedStructuralTopology& topology,
                                                         const SfamParameters& parameters,
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

std::vector<BondedTerm> SfamPotentialTermsGenerator::getBondedTerms() {
  return PotentialTermsHelper::getBondedTerms(topology_, parameters_, atomTypesHolder_);
}

std::vector<AngleTerm> SfamPotentialTermsGenerator::getAngleTerms() {
  return PotentialTermsHelper::getAngleTerms(topology_, parameters_, atomTypesHolder_);
}

std::vector<DihedralTerm> SfamPotentialTermsGenerator::getDihedralTerms() {
  std::vector<DihedralTerm> dihedralList;

  int numDihedralsMissing = 0;
  for (const auto& dihedral : topology_.getDihedralContainer()) {
    auto atomType1 = atomTypesHolder_.getAtomType(dihedral.atom1);
    auto atomType2 = atomTypesHolder_.getAtomType(dihedral.atom2);
    auto atomType3 = atomTypesHolder_.getAtomType(dihedral.atom3);
    auto atomType4 = atomTypesHolder_.getAtomType(dihedral.atom4);
    std::vector<Dihedral> m;
    try {
      // returns vector of dihedral parameters for given atom types, since it is
      // possible there are more than one in case of Fourier series
      m = parameters_.getMMDihedrals(atomType1, atomType2, atomType3, atomType4);
      DihedralType dihedralType(atomType1, atomType2, atomType3, atomType4);
      for (auto& d : m) {
        DihedralTerm term(dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4, d, dihedralType);
        dihedralList.push_back(term);
      }
    }
    catch (const MMDihedralParametersNotAvailableException& e) {
      log_.debug << e.what() << ".\n--> Dihedrals of this form will be ignored." << Core::Log::endl;
      numDihedralsMissing++;
    }
  }
  if (numDihedralsMissing > 0)
    log_.output << "\nThere are parameters missing for some of the dihedrals. These will be ignored." << Core::Log::nl
                << Core::Log::endl;

  return dihedralList;
}

std::vector<ImproperDihedralTerm> SfamPotentialTermsGenerator::getImproperDihedralTerms() {
  std::vector<ImproperDihedralTerm> improperDihedralList;

  for (const auto& improperDihedral : topology_.getImproperDihedralContainer()) {
    auto atomTypeC = atomTypesHolder_.getAtomType(improperDihedral.centralAtom);
    auto atomType2 = atomTypesHolder_.getAtomType(improperDihedral.atom2);
    auto atomType3 = atomTypesHolder_.getAtomType(improperDihedral.atom3);
    auto atomType4 = atomTypesHolder_.getAtomType(improperDihedral.atom4);
    ImproperDihedralType improperDihedralType(atomTypeC, atomType2, atomType3, atomType4);

    try {
      ImproperDihedral m = parameters_.getMMImproperDihedral(atomTypeC, atomType2, atomType3, atomType4);
      ImproperDihedralTerm term(improperDihedral.centralAtom, improperDihedral.atom2, improperDihedral.atom3,
                                improperDihedral.atom4, m, improperDihedralType);
      improperDihedralList.push_back(term);
    }
    catch (const MMImproperDihedralParametersNotAvailableException& e) {
      ImproperDihedral m;
      ImproperDihedralTerm term(improperDihedral.centralAtom, improperDihedral.atom2, improperDihedral.atom3,
                                improperDihedral.atom4, m, improperDihedralType);
      improperDihedralList.push_back(term);
    }
  }

  return improperDihedralList;
}

std::vector<DispersionTerm> SfamPotentialTermsGenerator::getDispersionTerms(bool applyCutoff) {
  std::vector<DispersionTerm> dispersionList;

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

      double factor = 1;
      if (type == -1)
        factor = scalingFactorForDispersionOneFourTerms_;

      auto atomTypeA = atomTypesHolder_.getAtomType(a);
      auto atomTypeB = atomTypesHolder_.getAtomType(b);
      auto c6 = static_cast<double>(parameters_.getC6(atomTypeA, atomTypeB));
      if (std::abs(c6) < 1e-6)
        throw std::runtime_error("Pair of atoms with indices " + std::to_string(a) + " and " + std::to_string(b) +
                                 " have a C6 coefficient of zero.");
      Dispersion disp(factor, c6);
      DispersionTerm term(a, b, disp, cutoff_);
      dispersionList.push_back(term);
    }
  }

  return dispersionList;
}

std::vector<RepulsionTerm> SfamPotentialTermsGenerator::getRepulsionTerms(bool applyCutoff) {
  std::vector<RepulsionTerm> repulsionList;

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

      double factor = 1;
      if (type == -1)
        factor = scalingFactorForRepulsionOneFourTerms_;
      Repulsion rep(factor);
      RepulsionTerm term(a, b, rep, cutoff_);
      repulsionList.push_back(term);
    }
  }

  return repulsionList;
}

std::vector<ElectrostaticTerm> SfamPotentialTermsGenerator::getElectrostaticTerms(bool applyCutoff) {
  // 0 = excluded, 1 = included, -1 = scaled
  Eigen::MatrixXi exclusionType = PotentialTermsHelper::getExclusionTypeMatrix(topology_, nAtoms_);
  return PotentialTermsHelper::getElectrostaticTerms(applyCutoff, cutoff_, scalingFactorForElectrostaticOneFourTerms_,
                                                     exclusionType, positions_);
}

std::vector<HydrogenBondTerm> SfamPotentialTermsGenerator::getHydrogenBondTerms() {
  std::vector<HydrogenBondTerm> hydrogenBondList;
  for (const auto& hydrogenBond : topology_.getHydrogenBondContainer()) {
    HydrogenBond hb;
    HydrogenBondTerm term(hydrogenBond.atom1, hydrogenBond.atom2, hydrogenBond.atom3, hb);
    hydrogenBondList.push_back(term);
  }

  return hydrogenBondList;
}

} // namespace MolecularMechanics
} // namespace Scine
