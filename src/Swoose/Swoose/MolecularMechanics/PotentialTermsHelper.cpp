/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "PotentialTermsHelper.h"
#include "AtomTypesHolder.h"
#include "MMExceptions.h"
#include "MMParameters.h"
#include "Topology/IndexedStructuralTopology.h"

namespace Scine {
namespace MolecularMechanics {
namespace PotentialTermsHelper {

Eigen::MatrixXi getExclusionTypeMatrix(const IndexedStructuralTopology& topology, int nAtoms) {
  // 0 = excluded, 1 = included, -1 = scaled
  Eigen::MatrixXi exclusionType = Eigen::MatrixXi::Ones(nAtoms, nAtoms);
  for (const auto& excluded : topology.getExcludedNonBondedContainer()) {
    exclusionType(excluded.atom1, excluded.atom2) = 0;
    exclusionType(excluded.atom2, excluded.atom1) = 0;
  }
  for (const auto& scaled : topology.getScaledNonBondedContainer()) {
    exclusionType(scaled.atom1, scaled.atom2) = -1;
    exclusionType(scaled.atom2, scaled.atom1) = -1;
  }

  return exclusionType;
}

std::vector<BondedTerm> getBondedTerms(const IndexedStructuralTopology& topology, const MMParameters& parameters,
                                       const AtomTypesHolder& atomTypesHolder) {
  std::vector<BondedTerm> bondList;

  for (const auto& bond : topology.getBondContainer()) {
    auto atomType1 = atomTypesHolder.getAtomType(bond.atom1);
    auto atomType2 = atomTypesHolder.getAtomType(bond.atom2);
    BondType bondType(atomType1, atomType2);
    try {
      Bond m = parameters.getMMBond(atomType1, atomType2);
      BondedTerm term(bond.atom1, bond.atom2, m, bondType);
      bondList.emplace_back(term);
    }
    catch (const MMBondParametersNotAvailableException& e) {
      Bond m; // bond is created empty, which will later throw the exception
      BondedTerm term(bond.atom1, bond.atom2, m, bondType);
      bondList.emplace_back(term);
    }
  }

  return bondList;
}

std::vector<AngleTerm> getAngleTerms(const IndexedStructuralTopology& topology, const MMParameters& parameters,
                                     const AtomTypesHolder& atomTypesHolder) {
  std::vector<AngleTerm> angleList;

  for (const auto& angle : topology.getAngleContainer()) {
    auto atomType1 = atomTypesHolder.getAtomType(angle.atom1);
    auto atomType2 = atomTypesHolder.getAtomType(angle.atom2);
    auto atomType3 = atomTypesHolder.getAtomType(angle.atom3);
    AngleType angleType(atomType1, atomType2, atomType3);
    try {
      Angle m = parameters.getMMAngle(atomType1, atomType2, atomType3);
      AngleTerm term(angle.atom1, angle.atom2, angle.atom3, m, angleType);
      angleList.push_back(term);
    }
    catch (const MMAngleParametersNotAvailableException& e) {
      Angle m;
      AngleTerm term(angle.atom1, angle.atom2, angle.atom3, m, angleType);
      angleList.push_back(term);
    }
  }

  return angleList;
}

std::vector<ElectrostaticTerm>
getElectrostaticTerms(bool applyCutoff, std::shared_ptr<double> cutoffRadius, double scalingFactorForOneFourTerms,
                      const Eigen::MatrixXi& exclusionTypeMatrix, const Utils::PositionCollection& positions) {
  std::vector<ElectrostaticTerm> electrostaticList;
  assert(positions.rows() == exclusionTypeMatrix.rows());
  assert(positions.rows() == exclusionTypeMatrix.cols());

  for (int a = 0; a < positions.rows(); ++a) {
    for (int b = 0; b < a; ++b) {
      if (applyCutoff) {
        const auto& R = (positions.row(b) - positions.row(a)).norm();
        if (R > *cutoffRadius)
          continue;
      }

      auto type = exclusionTypeMatrix(a, b);
      if (type == 0)
        continue;

      double factor = 1;
      if (type == -1)
        factor = scalingFactorForOneFourTerms;
      Electrostatic m(factor);
      ElectrostaticTerm term(a, b, m, cutoffRadius);
      electrostaticList.push_back(term);
    }
  }

  return electrostaticList;
}

} // namespace PotentialTermsHelper
} // namespace MolecularMechanics
} // namespace Scine
