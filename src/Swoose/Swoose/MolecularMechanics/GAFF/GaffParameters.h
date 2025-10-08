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
#include <memory>
#include <vector>

namespace Scine {

namespace Utils {
class AtomCollection;
} // namespace Utils

namespace MolecularMechanics {
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
  std::vector<Dihedral> getMMDihedrals(const std::string& t1, const std::string& t2, const std::string& t3,
                                       const std::string& t4) const;
  /** @brief Get ImproperDihedrals for the four atom types t1, t2, t3 and t4 */
  std::vector<Dihedral> getMMImproperDihedrals(const std::string& central, const std::string& t2, const std::string& t3,
                                               const std::string& t4) const;

  /** @brief Get Bond for the two atom types t1 and t2 */
  virtual Bond getMMBond(std::string t1, std::string t2) const override;
  /** @brief Get Angle for the three atom types t1, t2 and t3 */
  virtual Angle getMMAngle(std::string t1, std::string t2, std::string t3) const override;
  /**
   * @brief Get Lennard-Jones for two atom types t1 and t2.
   */
  const LennardJonesPairParameters& getMMLennardJones(const std::string& t1, const std::string& t2) const;

  // These functions add certain parameters of the MM model
  void addLennardJones(const std::string& atomType, LennardJonesParameters lennardJonesParameters);
  void addDihedral(const DihedralType& dihedralType, DihedralParameters dihedralParameters);
  void addImproperDihedral(const ImproperDihedralType& improperDihedralType, DihedralParameters improperDihedralParameters);
  void setType2AtomClass(const std::map<std::string, std::string>& map);
  /**
   * @return Return True if there were no parameter set.
   */
  bool empty() const;
  /**
   * @brief += operator to combine two parameter sets without any check for duplicated parameters.
   * @param other The other parameter set.
   * @return The combined parameters.
   */
  GaffParameters operator+=(const GaffParameters& other);
  /**
   * @brief Addition operator to combine two parameter sets.
   * @param lhs The left parameter set.
   * @param rhs The right parameter set.
   * @return The combined parameters.
   */
  friend GaffParameters operator+(GaffParameters lhs, const GaffParameters& rhs);
  /**
   * @brief Getter for the Lennard Jones parameters.
   * @return The atom type to LennardJonesParameter map.
   */
  const std::unordered_map<std::string, LennardJonesParameters>& getLennardJonesParameters() const;
  /**
   * @brief Setter for the 1-4 electrostatic scaling.
   * @param scaling The scaling factor.
   */
  void set14ElectrostaticScaling(double scaling);
  /**
   * @brief Setter for the 1-4 Lennard Jones scaling.
   * @param scaling The scaling factor.
   */
  void set14LennardJonesScaling(double scaling);
  /**
   * @brief Getter for the 1-4 electrostatic scaling.
   * @return The scaling factor.
   */
  double get14ElectrostaticScaling() const;
  /**
   * @brief Getter for the 1-4 Lennard Jones scaling.
   * @return The scaling factor.
   */
  double get14LennardJonesScaling() const;

 private:
  std::multimap<DihedralType, DihedralParameters> dihedrals_;
  std::multimap<ImproperDihedralType, DihedralParameters> improperDihedrals_;
  std::vector<std::vector<std::shared_ptr<LennardJonesPairParameters>>> lennardJonesPairs_;
  std::unordered_map<std::string, LennardJonesParameters> lennardJonesParameters_;
  std::map<std::string, std::string> type2AtomClass_;

  const std::string& getAtomClassFromType(const std::string& type) const;
  std::vector<Dihedral> resolveMMDihedrals(const std::string& t1, const std::string& t2, const std::string& t3,
                                           const std::string& t4) const;
  std::vector<Dihedral> resolveMMImproperDihedrals(const std::string& central, const std::string& t2,
                                                   const std::string& t3, const std::string& t4) const;

  double scalingFactorForElectrostaticOneFourTerms_ = 0.5;
  double scalingFactorForLennardJonesOneFourTerms_ = 0.5;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_GAFFPARAMETERS_H
