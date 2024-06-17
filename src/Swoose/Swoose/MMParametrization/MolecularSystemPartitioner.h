/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_MOLECULARSYSTEMPARTITIONER_H
#define MMPARAMETRIZATION_MOLECULARSYSTEMPARTITIONER_H

#include "ParametrizationUtils/ConstrainedAtomsIdentifier.h"
#include <Utils/Typenames.h>
#include <deque>
#include <memory>

namespace Scine {

namespace Core {
struct Log;
} // namespace Core

namespace Molassembler {
class InterpretResult;
} // namespace Molassembler

namespace Utils {
class Settings;
class Atom;
class AtomCollection;
} // namespace Utils

namespace MMParametrization {
struct ParametrizationData;
class FragmentAnalyzer;

/**
 * @class MolecularSystemPartitioner MolecularSystemPartitioner.h
 * @brief Divides the full system into subsystems.
 */
class MolecularSystemPartitioner {
 public:
  /**
   * @brief Constructor.
   */
  MolecularSystemPartitioner(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings, Core::Log& log);
  /**
   * @brief Divides the full system into subsystems and optionally saturates them with hydrogen atoms.
   *
   *        This function also creates the mapping of all relevant atom pairs to their respective subsystem.
   */
  void divideSystem();

 private:
  /*
   * @brief The function called when the number of atoms is smaller than the threshold given below.
   *
   *        This function does not perform a system division, just prepares the data_ object accordingly.
   */
  void prepareDataForOneSubsystem();
  /*
   * @brief Performs the system division if the number of atoms is larger than the threshold given below.
   */
  void divideSystemIntoSubsystems();
  /*
   * @brief The maximum number of atoms for which there will be no division into subsystems.
   */
  int numberAtomsThreshold_;
  /*
   * @brief The radius of the spheres defining the subsystems.
   */
  double subsystemRadius_;
  /*
   * @brief Threshold of which bond order to still consider as a bond.
   */
  double bondOrderThreshold_;
  // The data used within all MM parametrization classes
  ParametrizationData& data_;
  // The settings
  std::shared_ptr<Utils::Settings> settings_;
  // The class identifying which atoms have to be constrained during the geometry optimizations
  std::unique_ptr<ConstrainedAtomsIdentifier> constrainedAtomsIdentifier_;
  // Logger
  Core::Log& log_;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_MOLECULARSYSTEMPARTITIONER_H
