/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OptimizationSetup.h"
#include "ParametrizationData.h"
#include <Swoose/MolecularMechanics/Interactions/DihedralTerm.h>
#include <Swoose/MolecularMechanics/Interactions/ImproperDihedralTerm.h>
#include <Swoose/MolecularMechanics/Parameters/DispersionC6Parameters.h>
#include <Swoose/Utilities/OptionNames.h>
#include <Swoose/Utilities/SettingsNames.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Constants.h>

namespace Scine {
namespace MMParametrization {

constexpr double OptimizationSetup::initialBondForceConstant_;
constexpr double OptimizationSetup::initialAngleForceConstant_;
constexpr double OptimizationSetup::initialDihedralHalfBarrierHeight_;
constexpr double OptimizationSetup::initialImproperDihedralForceConstantForNonPlanarGroups_;
constexpr double OptimizationSetup::parameterValueInUninitializedState_;

OptimizationSetup::OptimizationSetup(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings)
  : data_(data), settings_(settings) {
  fragmentDataDistributor_ = std::make_unique<FragmentDataDistributor>(data);
}

void OptimizationSetup::generateInitialParameters() {
  setEquilibriumValues();
  setAtomicCharges();
  setNonCovalentParameters();
  setC6Parameters();
  setConstantDihedralParameters();
  setInitialGuessForForceConstants();
}

// This function extracts the equilibrium values for the internal coordinates from a reference structure
// and adds those values to the force field. If more than one coordinate is described by one set of atom types,
// a mean value is computed and used as the equilibrium value in the force field.
void OptimizationSetup::setEquilibriumValues() {
  // Bonds
  std::multimap<MolecularMechanics::BondType, double> bondMap;
  for (auto& bond : data_.topology.getBondContainer()) {
    auto pos = atomIndicesToPositions({bond.atom1, bond.atom2}, 2);
    auto pos1 = pos[0];
    auto pos2 = pos[1];
    auto r = (pos2 - pos1).norm() * Utils::Constants::angstrom_per_bohr;
    auto bondType =
        MolecularMechanics::BondType(data_.atomTypes.getAtomType(bond.atom1), data_.atomTypes.getAtomType(bond.atom2));
    bondMap.emplace(bondType, r);
  }

  // Iterate over the multimap which holds all bonds in order to add every bond type
  // with its mean value to the parameters.
  for (auto uniqueKeyIter = bondMap.begin(), end = bondMap.end(); uniqueKeyIter != end;
       uniqueKeyIter = bondMap.upper_bound(uniqueKeyIter->first)) {
    auto bondType = uniqueKeyIter->first;
    double meanValue = getMeanValueForEquilibriumValue(bondMap, bondType);
    data_.parameters.addBond(bondType, MolecularMechanics::BondParameters(parameterValueInUninitializedState_, meanValue));
  }

  // Angles
  std::multimap<MolecularMechanics::AngleType, double> angleMap;
  for (auto& angle : data_.topology.getAngleContainer()) {
    auto pos = atomIndicesToPositions({angle.atom2, angle.atom1, angle.atom3}, 1);
    // Note, that the order has to match the one above so that posX means atomX.
    auto pos1 = pos[1];
    auto pos2 = pos[0];
    auto pos3 = pos[2];

    Eigen::Vector3d a((pos1 - pos2));
    Eigen::Vector3d b((pos3 - pos2));

    double theta = acos((a.dot(b) / (a.norm() * b.norm()))) * Utils::Constants::degree_per_rad;
    auto angleType =
        MolecularMechanics::AngleType(data_.atomTypes.getAtomType(angle.atom1), data_.atomTypes.getAtomType(angle.atom2),
                                      data_.atomTypes.getAtomType(angle.atom3));
    angleMap.emplace(angleType, theta);
  }

  // Iterate over the multimap which holds all angles in order to add every angle type
  // with its mean value to the parameters.
  for (auto uniqueKeyIter = angleMap.begin(), end = angleMap.end(); uniqueKeyIter != end;
       uniqueKeyIter = angleMap.upper_bound(uniqueKeyIter->first)) {
    auto angleType = uniqueKeyIter->first;
    double meanValue = getMeanValueForEquilibriumValue(angleMap, angleType);
    data_.parameters.addAngle(angleType, MolecularMechanics::AngleParameters(parameterValueInUninitializedState_, meanValue));
  }

  // Improper dihedrals
  std::multimap<MolecularMechanics::ImproperDihedralType, double> improperDihedralMap;
  for (auto& improperDihedral : data_.topology.getImproperDihedralContainer()) {
    auto pos = atomIndicesToPositions(
        {improperDihedral.centralAtom, improperDihedral.atom2, improperDihedral.atom3, improperDihedral.atom4}, 4);
    auto pos1 = pos[0];
    auto pos2 = pos[1];
    auto pos3 = pos[2];
    auto pos4 = pos[3];

    Eigen::Vector3d F((pos1 - pos2));
    Eigen::Vector3d G((pos2 - pos3));
    Eigen::Vector3d H((pos4 - pos3));

    auto A = F.cross(G);
    auto B = H.cross(G);

    double theta = MolecularMechanics::ImproperDihedralTerm::getTheta(A, B, G) * Utils::Constants::degree_per_rad;

    auto improperDihedralType = MolecularMechanics::ImproperDihedralType(
        data_.atomTypes.getAtomType(improperDihedral.centralAtom), data_.atomTypes.getAtomType(improperDihedral.atom2),
        data_.atomTypes.getAtomType(improperDihedral.atom3), data_.atomTypes.getAtomType(improperDihedral.atom4));
    improperDihedralMap.emplace(improperDihedralType, theta);
  }

  // Iterate over the multimap which holds all improper dihedrals in order to add every improper dihedral type
  // with its mean value to the parameters.
  for (auto uniqueKeyIter = improperDihedralMap.begin(), end = improperDihedralMap.end(); uniqueKeyIter != end;
       uniqueKeyIter = improperDihedralMap.upper_bound(uniqueKeyIter->first)) {
    auto improperDihedralType = uniqueKeyIter->first;
    double meanValue = getMeanValueForEquilibriumValue(improperDihedralMap, improperDihedralType);

    // Set mean value of theta to 0 for planar groups
    // TODO: Set mean value of theta to 35.26 for tetrahedral groups? (commented out part)
    if (std::abs(meanValue) < MolecularMechanics::ImproperDihedral::getPlanarGroupThreshold()) {
      meanValue = 0.0;
    } //    else {
      //      meanValue = 35.26;
      //    }

    data_.parameters.addImproperDihedral(improperDihedralType, MolecularMechanics::ImproperDihedralParameters(
                                                                   parameterValueInUninitializedState_, meanValue));
  }
}

// This function adds the atomic charges to the force field.
void OptimizationSetup::setAtomicCharges() {
  assert(data_.atomicCharges.size() == data_.numberOfAtoms);
  assert(data_.atomTypes.size() == data_.numberOfAtoms);

  // Fill the charges multimap with all atomic charges along with the corresponding atom types
  std::multimap<std::string, double> chargesMap;
  for (int i = 0; i < data_.atomicCharges.size(); ++i) {
    chargesMap.emplace(data_.atomTypes.getAtomType(i), data_.atomicCharges.at(i));
  }

  // Iterate over the multimap which holds all atomic charges in order to add every atomic charge
  // with its mean value to the parameters.
  for (auto uniqueKeyIter = chargesMap.begin(), end = chargesMap.end(); uniqueKeyIter != end;
       uniqueKeyIter = chargesMap.upper_bound(uniqueKeyIter->first)) {
    auto atomType = uniqueKeyIter->first;
    double meanValue = getMeanValueForEquilibriumValue(chargesMap, atomType);
    data_.parameters.addCharge(atomType, meanValue);
  }
}

// This function adds the five non-covalent parameters (a1, s8, a2, beta, charges scaling factor) to the model.
void OptimizationSetup::setNonCovalentParameters() {
  std::vector<double> nonCovalentInitialParameters = {0.1, 4.6, 7.1, 7.4, 1.0};

  // Parameters were re-parametrized for use with the Turbomole Löwdin charges:
  bool turbomoleIsProgram = settings_->getString(SwooseUtilities::SettingsNames::referenceProgram) ==
                            SwooseUtilities::OptionNames::turbomoleOption;
  if (turbomoleIsProgram)
    nonCovalentInitialParameters = {0.2, 5.5, 6.2, 7.3, 0.9};

  data_.parameters.setNonCovalentParameters(nonCovalentInitialParameters);
}

// This function calculates the C6 parameters (as in D3) and adds them to the MM model.
void OptimizationSetup::setC6Parameters() {
  MolecularMechanics::DispersionC6Parameters::fillC6MatrixForCurrentStructure(data_.parameters, data_.fullStructure,
                                                                              data_.atomTypes);
}

// This function adds the periodicity and the phase shift to the dihedral part of the force field.
void OptimizationSetup::setConstantDihedralParameters() {
  for (const auto& dihedral : data_.topology.getDihedralContainer()) {
    int periodicity = getPeriodicity(dihedral.atom2, dihedral.atom3);
    double phaseShift = 0.0;
    double halfBarrierHeight = parameterValueInUninitializedState_;
    try {
      double ps = getPhaseShift(dihedral.atom2, dihedral.atom3, periodicity);
      phaseShift = ps;
    }
    catch (std::exception& e) {
      halfBarrierHeight = 0.0;
    }
    data_.parameters.addDihedral(MolecularMechanics::DihedralType("X", data_.atomTypes.getAtomType(dihedral.atom2),
                                                                  data_.atomTypes.getAtomType(dihedral.atom3), "X"),
                                 MolecularMechanics::DihedralParameters(halfBarrierHeight, phaseShift, periodicity));
  }
}

// This function adds the initial guess of the force constants to the force field.
void OptimizationSetup::setInitialGuessForForceConstants() {
  for (auto& bond : data_.parameters.getBonds()) {
    if (bond.second.getForceConstant() == parameterValueInUninitializedState_)
      bond.second.setForceConstant(initialBondForceConstant_);
  }

  for (auto& angle : data_.parameters.getAngles()) {
    if (angle.second.getForceConstant() == parameterValueInUninitializedState_)
      angle.second.setForceConstant(initialAngleForceConstant_);
  }

  for (auto& dihedral : data_.parameters.getDihedrals()) {
    if (dihedral.second.getHalfBarrierHeight() == parameterValueInUninitializedState_)
      dihedral.second.setHalfBarrierHeight(initialDihedralHalfBarrierHeight_);
  }

  for (auto& improperDihedral : data_.parameters.getImproperDihedrals()) {
    if (improperDihedral.second.getForceConstant() != parameterValueInUninitializedState_)
      continue;
    if (improperDihedral.second.getEquilibriumAngle() == 0.0) {
      improperDihedral.second.setForceConstant(getImproperDihedralForceConstantForPlanarGroups());
    }
    else {
      improperDihedral.second.setForceConstant(initialImproperDihedralForceConstantForNonPlanarGroups_);
    }
  }
}

// This function calculates the phase shift of the dihedral angles that describe the rotation
// around the atom1-atom2 bond.
// Implementation is based on the supporting information of
// the QuickFF parametrization scheme: J. Comput. Chem. 2015, 36, 1015-1027.
double OptimizationSetup::getPhaseShift(int atom1, int atom2, int periodicity) {
  auto fundamentalPeriod = (2 * M_PI / periodicity) * Utils::Constants::degree_per_rad;

  std::vector<double> images;
  enum interval { one, two, three };
  std::vector<interval> intervals;

  for (auto& dihedral : data_.topology.getDihedralContainer()) {
    if ((dihedral.atom2 == atom1 && dihedral.atom3 == atom2) || (dihedral.atom2 == atom2 && dihedral.atom3 == atom1)) {
      auto pos = atomIndicesToPositions({dihedral.atom2, dihedral.atom3, dihedral.atom1, dihedral.atom4}, 2);
      // Note, that the order has to match the one above so that posX means atomX.
      auto pos1 = pos[2];
      auto pos2 = pos[0];
      auto pos3 = pos[1];
      auto pos4 = pos[3];

      Eigen::Vector3d F((pos1 - pos2));
      Eigen::Vector3d G((pos2 - pos3));
      Eigen::Vector3d H((pos4 - pos3));

      auto A = F.cross(G);
      auto B = H.cross(G);

      double theta = std::abs(MolecularMechanics::DihedralTerm::getTheta(A, B, G) * Utils::Constants::degree_per_rad);
      while (theta > fundamentalPeriod) {
        theta -= fundamentalPeriod;
      }
      images.push_back(theta);
      if (theta <= fundamentalPeriod / 6.0) {
        intervals.push_back(one);
      }
      else if (theta <= 4.0 * fundamentalPeriod / 6.0) {
        intervals.push_back(two);
      }
      else {
        intervals.push_back(three);
      }
    }
  }

  if (std::all_of(intervals.begin() + 1, intervals.end(),
                  std::bind(std::equal_to<int>(), std::placeholders::_1, intervals.front()))) {
    double sumOfElements = 0.0;
    std::for_each(images.begin(), images.end(), [&](double n) { sumOfElements += n; });
    auto average = sumOfElements / images.size();
    average /= 10.0;
    average = round(average);
    average *= 10.0;
    return average * periodicity; // TODO: Is this correct? The paper says just take 0 if interval = one or three
  }
  // TODO: Improve the dihedrals which do not fall into the 'standard' category
  else {
    if (periodicity == 2) {
      // If interval three is in intervals, then the eq. angles are 0 and 180 degrees, otherwise 90 and -90 degrees.
      if (std::find(intervals.begin(), intervals.end(), three) != intervals.end())
        return 0.0;
      else
        return 90.0;
    }
    else {
      throw std::exception();
    }
  }
}

int OptimizationSetup::getPeriodicity(int atom1, int atom2) const {
  // Get number of neighbors of atom1 and atom2
  int n = data_.listsOfNeighbors[atom1].size();
  int m = data_.listsOfNeighbors[atom2].size();

  // Calculate and return periodicity
  return (n - 1) * (m - 1) / gcd(n - 1, m - 1);
}

int OptimizationSetup::gcd(int n1, int n2) const {
  return (n2 == 0) ? n1 : gcd(n2, n1 % n2);
}

std::vector<Eigen::RowVector3d> OptimizationSetup::atomIndicesToPositions(std::vector<int> indices,
                                                                          int numberOfInitialCandidateFragments) {
  assert(numberOfInitialCandidateFragments <= indices.size() &&
         "Number of specified initial candidate fragments is too large.");
  std::vector<Eigen::RowVector3d> vectorOfPositions;

  // First handle the case of no fragments, just one single molecular system.
  if (data_.vectorOfOptimizedStructures.size() == 1) {
    if (!data_.vectorOfOptimizedStructures[0])
      throw std::runtime_error("Equilibrium geometrical values could not be obtained.");
    for (const auto& i : indices)
      vectorOfPositions.emplace_back(data_.vectorOfOptimizedStructures[0]->at(i).getPosition());
    assert(vectorOfPositions.size() == indices.size());
  }
  // Now treat the case of more than one fragment.
  else {
    /*
     * Add the the first 'initialCandidateFragments' of the involved atoms to the vector of initial candidates
     */
    std::vector<int> initialCandidates;
    std::vector<int> candidates;
    for (int k = 0; k < numberOfInitialCandidateFragments; ++k) {
      initialCandidates.push_back(indices.at(k));
      candidates.push_back(indices.at(k));
    }

    for (const auto& index : initialCandidates) {
      fragmentDataDistributor_->updateCandidateFragments(index, candidates);
    }
    // We should have at least 4 candidates to ensure enough redundancy
    if (candidates.size() < 4) {
      for (const auto& index : initialCandidates) {
        fragmentDataDistributor_->updateCandidateFragmentsWithThirdShellNeighbors(index, candidates);
      }
    }

    // Loop over all of the atoms in the candidates vector and try to use their optimized structure.
    bool successful = false;
    for (const auto& i : candidates) {
      // The Hessian data should be obtained from the same fragment as these geometry-related data.
      if (!data_.vectorOfHessians.at(i))
        continue;
      if (!data_.vectorOfOptimizedStructures.at(i))
        throw std::runtime_error("No optimized structure available despite available Hessian. Fragment index: " +
                                 std::to_string(i));

      for (const auto& j : indices) {
        // Get index of j in the fragment i
        auto indexInFragment =
            std::distance(data_.atomIndexMapping[i].begin(),
                          std::find(data_.atomIndexMapping[i].begin(), data_.atomIndexMapping[i].end(), j));
        // Add position of j from the optimized structure i to the vector of positions
        try {
          vectorOfPositions.emplace_back(data_.vectorOfOptimizedStructures[i]->at(indexInFragment).getPosition());
        }
        catch (const std::out_of_range&) {
          throw std::runtime_error(
              "Error while obtaining equilibrium geometrical values from optimized fragments. A "
              "recognized bond/angle/dihedral is not present within the fragment it is supposed to be in.");
        }
      }
      successful = true;
      break;
    }
    if (!successful)
      throw std::runtime_error("Equilibrium geometrical values could not be obtained.");
    assert(vectorOfPositions.size() == indices.size());
  }

  return vectorOfPositions;
}

template<typename M, typename T>
double OptimizationSetup::getMeanValueForEquilibriumValue(const M& map, const T& parameterType) const {
  auto mapIteratorType = map.equal_range(parameterType);
  double numValues = std::distance(mapIteratorType.first, mapIteratorType.second);
  double totalValue = 0.0;

  for (auto it = mapIteratorType.first; it != mapIteratorType.second; it++)
    totalValue += it->second;

  return totalValue / numValues;
}

double OptimizationSetup::getImproperDihedralForceConstantForPlanarGroups() {
  return improperDihedralForceConstantForPlanarGroups_;
}

double OptimizationSetup::getInitialImproperDihedralForceConstantForNonPlanarGroups() {
  return initialImproperDihedralForceConstantForNonPlanarGroups_;
}

double OptimizationSetup::getInitialDihedralHalfBarrierHeight() {
  return initialDihedralHalfBarrierHeight_;
}

double OptimizationSetup::getInitialAngleForceConstant() {
  return initialAngleForceConstant_;
}

double OptimizationSetup::getInitialBondForceConstant() {
  return initialBondForceConstant_;
}

} // namespace MMParametrization
} // namespace Scine
