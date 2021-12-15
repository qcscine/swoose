/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMREGIONSELECTOR_H
#define SWOOSE_QMMM_QMREGIONSELECTOR_H

#include <Core/BaseClasses/ObjectWithLog.h>
#include <Utils/Geometry/AtomCollection.h>
#include <memory>
#include <vector>

namespace Scine {

namespace Utils {
class BondOrderCollection;
class Settings;
} // namespace Utils

namespace Qmmm {

/**
 * @struct QmmmModel QmRegionSelector.h
 * @brief Holds the information about one QM/MM model: its structure,
 *        its QM atom indices, and its charge and spin multiplicity.
 */
struct QmmmModel {
  Utils::AtomCollection structure;
  std::vector<int> qmAtomIndices;
  int molecularCharge;
  int spinMultiplicity;
};

/**
 * @class QmRegionSelector QmRegionSelector.h
 * @brief Automated selection of the QM region given a central atom around which it should be constructed
 *        and other settings.
 */
class QmRegionSelector : public Core::ObjectWithLog {
 public:
  /// @brief Constructor.
  QmRegionSelector();

  /**
   * @brief Generates the optimal QM region automatically.
   * @param fullSystem The molecular structure of the full system.
   */
  void generateQmRegion(const Utils::AtomCollection& fullSystem);

  /**
   * @brief Getter for the generated QM region as a Utils::AtomCollection.
   */
  Utils::AtomCollection getQmRegionStructure() const;

  /**
   * @brief Getter for the indices of the atoms in the generated QM region.
   */
  std::vector<int> getQmRegionIndices() const;

  /**
   * @brief Getter for the molecular charge and multiplicity of the generated QM region.
   */
  std::pair<int, int> getQmRegionChargeAndMultiplicity() const;

  /**
   * @brief Accessor for the settings.
   * @return Utils::Settings& The settings.
   */
  Utils::Settings& settings();

  /**
   * @brief Constant accessor for the settings.
   * @return const Utils::Settings& The settings.
   */
  const Utils::Settings& settings() const;

 private:
  // Helper function that reads in the system connectivity from the connectivity file specified in the settings.
  Utils::BondOrderCollection getBondOrders(const Utils::AtomCollection& structure) const;
  // The settings.
  std::unique_ptr<Utils::Settings> settings_;
  // The QM/MM model candidates.
  std::vector<QmmmModel> qmmmModelCandidates_;
  // The QM/MM reference models.
  std::vector<QmmmModel> qmmmReferenceModels_;
  // The selected QM region. Initialized as -1, which is the state that no selection has been made yet.
  int selectedQmRegionIndex_ = -1;
};

class QmRegionHasNotBeenSelectedException : public std::runtime_error {
 public:
  explicit QmRegionHasNotBeenSelectedException() : std::runtime_error("The QM Region has not been generated yet.") {
  }
};

} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMREGIONSELECTOR_H
