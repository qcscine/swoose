/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_INTERACTIONTERMBASE_H
#define MOLECULARMECHANICS_INTERACTIONTERMBASE_H

namespace Scine {
namespace MolecularMechanics {

/**
 * @class InteractionTermBase InteractionTermBase.h
 * @brief Base class for all interaction terms.
 */
class InteractionTermBase {
 public:
  /**
   * @brief Disable this term. If it is disabled, the interaction is included in the MM model.
   *        Needed for the QM/MM calculator to switch specific interactions off and on.
   */
  void disable() {
    disabled_ = true;
  };

  /**
   * @brief Enable this term. If it is enabled, the interaction is included in the MM model.
   *        Needed for the QM/MM calculator to switch specific interactions off and on.
   */
  void enable() {
    disabled_ = false;
  };

 protected:
  bool disabled_ = false;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_INTERACTIONTERMBASE_H
