/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "IndexedStructuralTopology.h"

namespace Scine {
namespace MolecularMechanics {

void IndexedStructuralTopology::addBond(int a1, int a2) {
  bonds_.emplace_back(IndexedStructuralBond(a1, a2));
}

void IndexedStructuralTopology::addAngle(int a1, int a2, int a3) {
  angles_.emplace_back(IndexedStructuralAngle(a1, a2, a3));
}

void IndexedStructuralTopology::addDihedral(int a1, int a2, int a3, int a4) {
  dihedrals_.emplace_back(IndexedStructuralDihedral(a1, a2, a3, a4));
}

void IndexedStructuralTopology::addImproperDihedral(int central, int a2, int a3, int a4) {
  improperDihedrals_.emplace_back(IndexedStructuralImproperDihedral(central, a2, a3, a4));
}

void IndexedStructuralTopology::addExcludedNonBonded(int a1, int a2) {
  excludedNB_.emplace_back(StructuralExcludedNonBonded(a1, a2));
}

void IndexedStructuralTopology::addScaledNonBonded(int a1, int a2) {
  scaledNB_.emplace_back(IndexedStructuralScaledNonBonded(a1, a2));
}

void IndexedStructuralTopology::addHydrogenBond(int a1, int a2, int a3) {
  hydrogenBonds_.emplace_back(IndexedStructuralHydrogenBond(a1, a2, a3));
}

const std::vector<IndexedStructuralBond>& IndexedStructuralTopology::getBondContainer() const {
  return bonds_;
}

std::vector<IndexedStructuralBond>& IndexedStructuralTopology::getBondContainer() {
  return bonds_;
}

const std::vector<IndexedStructuralAngle>& IndexedStructuralTopology::getAngleContainer() const {
  return angles_;
}

std::vector<IndexedStructuralAngle>& IndexedStructuralTopology::getAngleContainer() {
  return angles_;
}

const std::vector<IndexedStructuralDihedral>& IndexedStructuralTopology::getDihedralContainer() const {
  return dihedrals_;
}

std::vector<IndexedStructuralDihedral>& IndexedStructuralTopology::getDihedralContainer() {
  return dihedrals_;
}

const std::vector<IndexedStructuralImproperDihedral>& IndexedStructuralTopology::getImproperDihedralContainer() const {
  return improperDihedrals_;
}

std::vector<IndexedStructuralImproperDihedral>& IndexedStructuralTopology::getImproperDihedralContainer() {
  return improperDihedrals_;
}

const std::vector<StructuralExcludedNonBonded>& IndexedStructuralTopology::getExcludedNonBondedContainer() const {
  return excludedNB_;
}

const std::vector<IndexedStructuralScaledNonBonded>& IndexedStructuralTopology::getScaledNonBondedContainer() const {
  return scaledNB_;
}

const std::vector<IndexedStructuralHydrogenBond>& IndexedStructuralTopology::getHydrogenBondContainer() const {
  return hydrogenBonds_;
}

void IndexedStructuralTopology::clearHydrogenBonds() {
  hydrogenBonds_.clear();
}

} // namespace MolecularMechanics
} // namespace Scine
