/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_DATAMANAGER_H
#define MMPARAMETRIZATION_DATAMANAGER_H

#include <Swoose/MolecularMechanics/AtomTypesHolder.h>
#include <Swoose/MolecularMechanics/SFAM/SfamParameters.h>
#include <Swoose/MolecularMechanics/Topology/IndexedStructuralTopology.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Typenames.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <list>
#include <memory>

namespace Scine {
namespace MMParametrization {
/**
 * @struct ParametrizationData ParametrizationData.h
 * @brief This struct holds all objects used inside the MM parametrization algorithm.
 */
struct ParametrizationData {
  /**
   * @brief The structure of the entire system to parametrize.
   */
  Utils::AtomCollection fullStructure;
  /**
   * @brief Number of atoms in the full structure.
   */
  int numberOfAtoms{0};
  /**
   * @brief Full Hessian matrix of the system that only contains the subblocks needed for the parametrization procedure.
   */
  Eigen::SparseMatrix<double> fullHessian;
  /**
   * @brief Vector of unique pointers to dense Eigen matrices that represent the Hessians
   *        of the subsystems of the parametrization.
   */
  std::vector<std::unique_ptr<Utils::HessianMatrix>> vectorOfHessians;
  /**
   * @brief Vector of unique pointers to molecular structures
   *        that represent the subsystems of the parametrization.
   */
  std::vector<std::unique_ptr<Utils::AtomCollection>> vectorOfStructures;
  /**
   * @brief Vector of unique pointers to molecular structures
   *        that represent the optimized subsystems of the parametrization.
   */
  std::vector<std::unique_ptr<Utils::AtomCollection>> vectorOfOptimizedStructures;
  /**
   * @brief Vector of charge and spin multiplicity pairs for each molecular structure
   *        representing a subsystem of the parametrization.
   */
  std::vector<std::pair<int, int>> vectorOfChargesAndMultiplicities;
  /**
   * @brief Vector of unique pointers to bond order matrices for for each molecular structure
   *        representing a subsystem of the parametrization.
   */
  std::vector<std::unique_ptr<Utils::BondOrderCollection>> vectorOfBondOrderCollections;
  /**
   * @brief Atomic charges for the whole system calculated from fragments.
   */
  std::vector<double> atomicCharges;
  /**
   * @brief Vector that holds the atomic charges for each fragment of the system.
   */
  std::vector<std::vector<double>> atomicChargesForEachFragment;
  /**
   * @brief The connectivity of the system. It is a vector of a list of neighbor atom indices for each atom.
   */
  std::vector<std::list<int>> listsOfNeighbors;
  /**
   * @brief Topology of the system.
   */
  MolecularMechanics::IndexedStructuralTopology topology;
  /**
   * @brief The SFAM molecular mechanics parameters.
   */
  MolecularMechanics::SfamParameters parameters;
  /**
   * @brief The atom types.
   */
  MolecularMechanics::AtomTypesHolder atomTypes;
  /**
   * @brief Vector that holds the pairs of atom indices corresponding to the partial Hessian blocks
   *        relevant for the current parameter's optimization.
   */
  std::vector<std::pair<int, int>> vectorOfHessianSubmatrixIndices;
  /**
   * @brief Vector that holds the atom types of the parameter currently optimized.
   */
  std::vector<std::string> vectorOfAtomTypesForParameter;
  /**
   * @brief A bond order matrix of the full system.
   */
  Utils::BondOrderCollection bondOrders;
  /**
   * @brief A map containing the indices of atoms and their formal charge in the full system.
   */
  std::map<int, int> formalCharges;
  /**
   * @brief A map containing the indices of atoms and the number of unpaired electrons that can
   *        be assigned to that atom.
   */
  std::map<int, int> unpairedElectrons;
  /**
   * @brief Vector of all fragments containing a vector of indices that correspond to the indices
   *        of the atoms inside the fragment in the full system.
   */
  std::vector<std::vector<int>> atomIndexMapping;
  /**
   * @brief Contains the indices of the atoms to be constrained during a geometry optimization
   *        for each fragment.
   */
  std::vector<std::vector<int>> constrainedAtoms;
  /**
   * @brief A vector of fragment indices which are superfluous and for which the reference data
   *        does not need to be calculated.
   */
  std::vector<int> superfluousFragments;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_DATAMANAGER_H
