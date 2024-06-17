/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_GAFFPARAMETERS_H
#define MOLECULARMECHANICS_GAFFPARAMETERS_H

#include "../MMParameters.h"
#include "../Parameters/AngleParameters.h"
#include "../Parameters/BondParameters.h"
#include "../Parameters/DihedralParameters.h"
#include "../Parameters/ImproperDihedralParameters.h"
#include "../Parameters/LennardJonesParameters.h"
#include "../Topology/AngleType.h"
#include "../Topology/BondType.h"
#include "../Topology/DihedralType.h"
#include "../Topology/ImproperDihedralType.h"
#include <Eigen/Core>
#include <map>
#include <vector>

namespace Scine {

namespace Utils {
class AtomCollection;
} // namespace Utils

namespace MolecularMechanics {
class AtomTypesHolder;
/**
 * @class GaffParameters GaffParameters.h
 * @brief Class containing the parameters for the GAFF model obtained after parsing a GAFF parameter file.
 *        The angle and bond parameters are handled by the base class MMParameters.
 */
class GaffParameters : public MMParameters {
 public:
  /**
   * @brief Get Dihedrals for the four atom types t1, t2, t3 and t4.
   *        (in principle there could be more than one (Fourier series))
   */
  std::vector<Dihedral> getMMDihedrals(std::string t1, std::string t2, std::string t3, std::string t4) const;
  /** @brief Get ImproperDihedrals for the four atom types t1, t2, t3 and t4 */
  std::vector<Dihedral> getMMImproperDihedrals(std::string central, std::string t2, std::string t3, std::string t4) const;
  /**
   * @brief Get Lennard-Jones for two atom types t1 and t2.
   */
  LennardJones getMMLennardJones(std::string t1, std::string t2, double scalingFactor) const;

  // These functions add certain parameters of the MM model
  void addLennardJones(std::string atomType, LennardJonesParameters lennardJonesParameters);
  void addDihedral(DihedralType dihedralType, DihedralParameters dihedralParameters);
  void addImproperDihedral(ImproperDihedralType improperDihedralType, DihedralParameters improperDihedralParameters);
  /**
   * @return Return True if there were no parameter set.
   */
  bool empty();

 private:
  std::multimap<DihedralType, DihedralParameters> dihedrals_;
  std::multimap<ImproperDihedralType, DihedralParameters> improperDihedrals_;
  std::map<std::string, LennardJonesParameters> lennardJonesPairs_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_GAFFPARAMETERS_H
