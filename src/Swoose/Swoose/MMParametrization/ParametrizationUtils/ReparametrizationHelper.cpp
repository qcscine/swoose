/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ReparametrizationHelper.h"
#include "../ParametrizationData.h"
#include <Core/Log.h>
#include <Swoose/MolecularMechanics/SFAM/SfamParameterParser.h>
#include <Swoose/MolecularMechanics/Topology/IndexedStructuralTopology.h>

namespace Scine {
namespace MMParametrization {

ReparametrizationHelper::ReparametrizationHelper(ParametrizationData& data, Core::Log& log) : data_(data), log_(log) {
}

void ReparametrizationHelper::parseProvidedParameters(const std::string& parameterFile) {
  log_.output << "Parsing provided parameter file '" << parameterFile << "' to extract existing parameters..."
              << Core::Log::endl;
  MolecularMechanics::SfamParameterParser parser(parameterFile, data_.atomTypes);
  data_.parameters = *parser.parseParameters();
}

void ReparametrizationHelper::manipulateTopology() {
  // Next two lines are to make the statements below shorter:
  using namespace MolecularMechanics;
  auto atomTypes = [&](int i) { return data_.atomTypes.getAtomType(i); };

  // Bonds:
  auto& bonds = data_.topology.getBondContainer();
  auto& bondParams = data_.parameters.getBonds();
  bonds.erase(std::remove_if(bonds.begin(), bonds.end(),
                             [&](const IndexedStructuralBond& b) {
                               BondType t(atomTypes(b.atom1), atomTypes(b.atom2));
                               bool found = bondParams.find(t) != bondParams.end();
                               if (!found)
                                 relevantFragments_.insert(b.atom1); // TODO: add more?
                               return found;
                             }),
              bonds.end());

  // Angles:
  auto& angles = data_.topology.getAngleContainer();
  auto& angleParams = data_.parameters.getAngles();
  angles.erase(std::remove_if(angles.begin(), angles.end(),
                              [&](const IndexedStructuralAngle& a) {
                                AngleType t(atomTypes(a.atom1), atomTypes(a.atom2), atomTypes(a.atom3));
                                bool found = angleParams.find(t) != angleParams.end();
                                if (!found)
                                  relevantFragments_.insert(a.atom2); // TODO: add more?
                                return found;
                              }),
               angles.end());

  // Dihedrals:
  auto& dihedrals = data_.topology.getDihedralContainer();
  auto& dihedralParams = data_.parameters.getDihedrals();
  dihedrals.erase(std::remove_if(dihedrals.begin(), dihedrals.end(),
                                 [&](const IndexedStructuralDihedral& d) {
                                   DihedralType t("X", atomTypes(d.atom2), atomTypes(d.atom3), "X");
                                   bool found = dihedralParams.find(t) != dihedralParams.end();
                                   if (!found)
                                     relevantFragments_.insert(d.atom2); // TODO: add more?
                                   return found;
                                 }),
                  dihedrals.end());

  // Improper dihedrals:
  auto& impropers = data_.topology.getImproperDihedralContainer();
  auto& improperParams = data_.parameters.getImproperDihedrals();
  impropers.erase(std::remove_if(impropers.begin(), impropers.end(),
                                 [&](const IndexedStructuralImproperDihedral& impr) {
                                   ImproperDihedralType t(atomTypes(impr.centralAtom), atomTypes(impr.atom2),
                                                          atomTypes(impr.atom3), atomTypes(impr.atom4));
                                   bool found = improperParams.find(t) != improperParams.end();
                                   if (!found)
                                     relevantFragments_.insert(impr.centralAtom); // TODO: add more?
                                   return found;
                                 }),
                  impropers.end());

  // Atomic charges:
  auto& atomicChargesParams = data_.parameters.getCharges();
  for (int i = 0; i < data_.numberOfAtoms; ++i) {
    bool found = atomicChargesParams.find(atomTypes(i)) != atomicChargesParams.end();
    if (!found)
      relevantFragments_.insert(i);
  }
}

bool ReparametrizationHelper::isRelevantFragment(int fragmentIndex) const {
  return relevantFragments_.find(fragmentIndex) != relevantFragments_.end();
}

} // namespace MMParametrization
} // namespace Scine
