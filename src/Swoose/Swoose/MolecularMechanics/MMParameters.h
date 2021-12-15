/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_MMPARAMETERS_H
#define MOLECULARMECHANICS_MMPARAMETERS_H

#include "Parameters/AngleParameters.h"
#include "Parameters/BondParameters.h"
#include "Topology/AngleType.h"
#include "Topology/BondType.h"
#include <map>

namespace Scine {

namespace Utils {
class AtomCollection;
} // namespace Utils

namespace MolecularMechanics {
/**
 * @class MMParameters MMParameters.h
 * @brief Base class for the MM parameters classes, which holds and manages the bond and angle parameters as these
 *        are identically handled for SFAM and GAFF. Functions are still labeled as virtual such that these can be
 *        overwritten.
 */
class MMParameters {
 public:
  /** @brief Get Bond for the two atom types t1 and t2 */
  virtual Bond getMMBond(std::string t1, std::string t2) const;
  /** @brief Get Angle for the three atom types t1, t2 and t3 */
  virtual Angle getMMAngle(std::string t1, std::string t2, std::string t3) const;

  /** @brief Adds a bond to the bonds_ container */
  virtual void addBond(BondType bondType, BondParameters bondParameters);
  /** @brief Adds an angle to the angles_ container */
  virtual void addAngle(AngleType angleType, AngleParameters angleParameters);

 protected:
  std::map<BondType, BondParameters> bonds_;
  std::map<AngleType, AngleParameters> angles_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_MMPARAMETERS_H