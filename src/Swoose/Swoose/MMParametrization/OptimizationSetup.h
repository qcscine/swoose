/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_OPTIMIZATIONSETUP_H
#define MMPARAMETRIZATION_OPTIMIZATIONSETUP_H

#include "ParametrizationUtils/FragmentDataDistributor.h"
#include <Eigen/Core>
#include <memory>
#include <vector>

namespace Scine {

namespace Utils {
class Settings;
} // namespace Utils

namespace MMParametrization {
struct ParametrizationData;

/**
 * @class OptimizationSetup OptimizationSetup.h
 * @brief This class sets up the initial MM parameters. The force constants will be optimized after that.
 */
class OptimizationSetup {
 public:
  /**
   * @brief Constructor.
   */
  OptimizationSetup(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings);
  /**
   * @brief This function generates the initial MM parameters.
   */
  void generateInitialParameters();
  /**
   * @brief Getter for the constant value that is used for the force constants of improper dihedrals.
   */
  static double getImproperDihedralForceConstantForPlanarGroups();
  /**
   * @brief Getter for the initial value that is used for the force constants of improper dihedrals of non-planar groups.
   */
  static double getInitialImproperDihedralForceConstantForNonPlanarGroups();
  /**
   * @brief Getter for the initial value that is used for the force constants of bonds.
   */
  static double getInitialBondForceConstant();
  /**
   * @brief Getter for the initial value that is used for the force constants of angles.
   */
  static double getInitialAngleForceConstant();
  /**
   * @brief Getter for the initial value that is used for the force constants of dihedral angles.
   */
  static double getInitialDihedralHalfBarrierHeight();

 private:
  /*
   * @brief Converts a vector of atom indices to a vector of the positions of these atoms in
   *        a reasonable optimized fragment. The second argument is the number of indices from that
   *        given vector (starting from front) which should be used as initial candidate indices to build
   *        the vector of candidate fragments from which to take the reference data. For example, if the indices of
   *        an angle are {atom1, atom2, atom3} and numberOfInitialCandidateFragments is 2, atom1 and atom2 build
   *        the initial candidate fragments vector which is then updated by the FragmentDataDistributor class
   *        to include their neighbors.
   */
  std::vector<Eigen::RowVector3d> atomIndicesToPositions(std::vector<int> indices, int numberOfInitialCandidateFragments);
  /*
   * @brief Extracts all equilibrium values in the MM parameters from the reference structure.
   */
  void setEquilibriumValues();
  /*
   * @brief Adds atomic charges to the MM parameters.
   */
  void setAtomicCharges();
  /*
   * @brief Adds non-covalent parameters to the MM parameters.
   */
  void setNonCovalentParameters();
  /*
   * @brief Adds C6 parameters to the MM parameters.
   */
  void setC6Parameters();
  /*
   * @brief Adds periodicity and phase shift for dihedral angles to MM parameters.
   */
  void setConstantDihedralParameters();
  /*
   * @brief Creates the initial guess for the force constants and barrier heights.
   */
  void setInitialGuessForForceConstants();
  /*
   * @brief This function evaluates the mean value for a given parameter given a map containing the parameter types
   *        and their values as well as a current parameter type to consider.
   */
  template<typename M, typename T>
  double getMeanValueForEquilibriumValue(const M& map, const T& parameterType) const;
  // This function calculates the phase shift of a dihedral angle.
  double getPhaseShift(int atom1, int atom2, int periodicity);
  // This function calculates the periodicity of a dihedral angle.
  int getPeriodicity(int atom1, int atom2) const;
  // Greatest common denominator function.
  int gcd(int n1, int n2) const;
  // Constant value that is used for the force constants of improper dihedrals.
  static constexpr double improperDihedralForceConstantForPlanarGroups_ = 50.0; // TODO: What is a good value?
  // The data used within all MM parametrization classes
  ParametrizationData& data_;
  // The settings
  std::shared_ptr<Utils::Settings> settings_;
  // Pointer to an instance of the fragment data distributor class needed to find fragment candidates to get data from.
  std::unique_ptr<FragmentDataDistributor> fragmentDataDistributor_;
  // Initial force constant values
  static constexpr double initialBondForceConstant_ = 700.0;       // Unit: kcal/(mol*A^2)
  static constexpr double initialAngleForceConstant_ = 120.0;      // Unit: kcal/(mol*rad^2)
  static constexpr double initialDihedralHalfBarrierHeight_ = 0.5; // Unit kcal/mol
  // TODO: Does it make sense to start with 0.0?
  static constexpr double initialImproperDihedralForceConstantForNonPlanarGroups_ = 0.0; // Unit: kcal/(mol*rad^2)
  /*
   * When a parameter is created, it does not get assigned an initial value
   * for the force constant/half barrier height right away, but only in
   * a later function. Therefore, this constant is introduced to keep
   * track whether a correct initial value has been set already for
   * a given parameter. This is important, because if a parameter already
   * exists, we do not want to overwrite its value.
   */
  static constexpr double parameterValueInUninitializedState_ = -1.0;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_OPTIMIZATIONSETUP_H
