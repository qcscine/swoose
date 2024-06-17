/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SubsystemGenerator.h"
#include "FragmentAnalyzer.h"
#include "FragmentationHelper.h"
#include <Core/Log.h>
#include <Molassembler/Graph.h>
#include <Molassembler/Interpret.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondDetector.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>

namespace Scine {
namespace SwooseUtilities {

namespace {
static constexpr const char* failedAttemptStructureSmallestRadiusFilename = "failed_attempt_small.xyz";
static constexpr const char* failedAttemptStructureLargestRadiusFilename = "failed_attempt_large.xyz";
} // namespace

SubsystemGenerator::SubsystemGenerator(const Utils::AtomCollection& fullStructure,
                                       const Utils::BondOrderCollection& bondOrders, FragmentAnalyzer& fragmentAnalyzer,
                                       double bondOrderThreshold, int maximumSubsystemSize, Core::Log& log,
                                       int randomSeed, double probabilityToDivideBond)
  : fullStructure_(fullStructure),
    bondOrders_(bondOrders),
    fragmentAnalyzer_(fragmentAnalyzer),
    bondOrderThreshold_(bondOrderThreshold),
    maximumSubsystemSize_(maximumSubsystemSize),
    randomEngine_(std::make_shared<std::mt19937>(randomSeed)),
    log_(log),
    probabilityToDivideBond_(probabilityToDivideBond) {
  // Fill up the vector containing the subgraph sizes, which is needed later.
  subgraphSizes_ = FragmentationHelper::calculateSubgraphSizes(fullStructure_, bondOrders_);
  // Generate lists of neighbors from bond order matrix.
  listsOfNeighbors_ = SwooseUtilities::TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix(
      fullStructure_.size(), bondOrders_, bondOrderThreshold_);
}

Utils::AtomCollection SubsystemGenerator::generateSubsystem(int centralAtomIndex, std::vector<int>& atomIndexMapping,
                                                            double subsystemRadius) {
  // Stuff to initiate before while loop
  double additionToRadius = 0.0;
  Utils::AtomCollection subsystem;
  Utils::AtomCollection firstAttemptSubsystem;
  Utils::Atom centralAtom = fullStructure_.at(centralAtomIndex);
  bool unsuccessful = true;
  int attemptCounter = 0;
  // while loop to increase subsystem radius as long as the system is valid
  while (unsuccessful) {
    attemptCounter++;
    if (attemptCounter > maxNumOfAttemptsForOneSubsystem_) {
      Utils::ChemicalFileHandler::write(failedAttemptStructureSmallestRadiusFilename, firstAttemptSubsystem);
      Utils::ChemicalFileHandler::write(failedAttemptStructureLargestRadiusFilename, subsystem);
      log_.output
          << Core::Log::nl
          << "The fragment around atom " + std::to_string(centralAtomIndex) + " could not be generated." << Core::Log::nl
          << Core::Log::nl << "For debugging, the generated structures of this fragment with the smallest and the largest initial radius can be found in \""
          << failedAttemptStructureSmallestRadiusFilename << "\" and \"" << failedAttemptStructureLargestRadiusFilename
          << "\"." << Core::Log::nl << Core::Log::nl << "Please check if you provided the correct atomic info file. "
          << Core::Log::endl;
      throw std::runtime_error("Error while trying to generate fragment with index " + std::to_string(centralAtomIndex) + ".");
    }
    tryGeneratingSensibleSubsystem(subsystem, centralAtom, centralAtomIndex, atomIndexMapping, subsystemRadius,
                                   additionToRadius, unsuccessful);
    if (attemptCounter == 1)
      firstAttemptSubsystem = subsystem;
  }
  return subsystem;
}

void SubsystemGenerator::tryGeneratingSensibleSubsystem(Utils::AtomCollection& subsystem, const Utils::Atom& centralAtom,
                                                        int atomIndex, std::vector<int>& atomIndexMapping,
                                                        const double& subsystemRadius, double& additionToRadius,
                                                        bool& unsuccessful) {
  // First generate a rough guess for the subsystem by strictly applying the cutoff radius
  Utils::AtomCollection preliminarySubsystem;
  std::vector<std::pair<Utils::Atom, int>> potentialBondingPartners;
  std::vector<int> indexMap = {atomIndex};
  preliminarySubsystem.push_back(centralAtom); // The first atom is the one that the subsystem is built around
  constexpr double devThresh = 2 * Utils::Constants::bohr_per_angstrom;
  int otherIndex = 0;
  for (const auto& otherAtom : fullStructure_) {
    if (atomIndex == otherIndex) {
      otherIndex++;
      continue;
    }
    // Calculate the distance between the two atoms
    auto distance = (otherAtom.getPosition() - centralAtom.getPosition()).norm();
    auto deviation = distance - (subsystemRadius + additionToRadius);
    // Add all atoms that are close enough, but not the same atom
    if (deviation < 0) {
      preliminarySubsystem.push_back(otherAtom);
      indexMap.push_back(otherIndex);
    }
    else if (deviation < devThresh) {
      potentialBondingPartners.emplace_back(std::make_pair(otherAtom, otherIndex));
    }
    otherIndex++;
  }

  subsystem.clear();
  subsystem.push_back(centralAtom); // The first atom is the one that the subsystem is built around
  auto bondOrders = Utils::BondDetector::detectBonds(preliminarySubsystem);

  auto result = Molassembler::Interpret::graphs(preliminarySubsystem, bondOrders);

  assert(preliminarySubsystem.size() == result.componentMap.size());

  Utils::AtomCollection atomsToAdd;   // In this AtomCollection the atoms to add to the subsystem are collected.
  std::vector<int> atomsToAddIndices; // This is a vector of atom indices of atoms that are already in the full
  // system and that are added to the subsystem later.
  std::deque<bool> isSaturatingAtom; // This is a deque that lists for all atoms in atomsToAdd whether it is an
  // atom placed in the subsystem for valence saturation.
  for (int i = 0; i < preliminarySubsystem.size(); ++i) {
    auto nAtoms = result.graphs[result.componentMap.apply(i).component].V();
    auto distance = (preliminarySubsystem.getPosition(i) - centralAtom.getPosition()).norm();

    // Should this atom be added to the subsystem?
    bool validAtom = nAtoms > 4 || subgraphSizes_.at(indexMap[i]) <= 4;

    if (validAtom) {
      if (distance > 1e-2) { // If it is the central atom, then it was already added.
        subsystem.push_back(preliminarySubsystem.at(i));
        for (const auto& pbp : potentialBondingPartners) {
          if (bondOrders_.getOrder(indexMap[i], pbp.second) > bondOrderThreshold_) {
            FragmentationHelper::addAtomsUpToReasonableCut(atomsToAdd, atomsToAddIndices, isSaturatingAtom, indexMap[i],
                                                           pbp.second, fullStructure_, listsOfNeighbors_,
                                                           probabilityToDivideBond_, randomEngine_);
          }
        }
      }
    }
    else {
      if (distance < 1e-2) { // Not a valid atom, but it is the central atom? Error!
        throw std::runtime_error("The fragment for atom " + std::to_string(atomIndex) + " could not be generated.");
      }
    }
  }

  // Add the atoms of the atomsToAdd structure to the final version of the subsystem -> First, non-saturating atoms
  FragmentationHelper::addMoreAtomsToSubsystem(atomsToAdd, subsystem, isSaturatingAtom, false);

  // The same for the saturating atoms
  FragmentationHelper::addMoreAtomsToSubsystem(atomsToAdd, subsystem, isSaturatingAtom, true);

  // Update index mapping information now since it is required in the analyzeFragment function
  FragmentationHelper::updateInformationForIndexMapping(subsystem, fullStructure_, atomIndexMapping);

  // Check whether the fragment is closed shell, if it is not, repeat procedure with a larger sphere
  bool valid = fragmentAnalyzer_.analyzeFragment(subsystem, atomIndexMapping);
  if (valid)
    unsuccessful = false;
  else {
    additionToRadius += 0.1; // add 0.1 bohr
  }

  // Check whether subsystem is too small or large
  FragmentationHelper::checkSizeOfSubsystem(subsystem.size(), additionToRadius, unsuccessful, atomIndex,
                                            fullStructure_.size(), maximumSubsystemSize_, log_);
}

} // namespace SwooseUtilities
} // namespace Scine
