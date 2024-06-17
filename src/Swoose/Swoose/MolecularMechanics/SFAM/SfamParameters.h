/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_SFAMPARAMETERS_H
#define MOLECULARMECHANICS_SFAMPARAMETERS_H

#include "../MMParameters.h"
#include "../Parameters/AngleParameters.h"
#include "../Parameters/BondParameters.h"
#include "../Parameters/DihedralParameters.h"
#include "../Parameters/ImproperDihedralParameters.h"
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
 * @class SfamParameters SfamParameters.h
 * @brief Class containing the parameters for SFAM's MM model obtained after parsing a SFAM parameter file.
 *        The angle and bond parameters are handled by the base class MMParameters.
 */
class SfamParameters : public MMParameters {
 public:
  /** @brief Checks whether the SFAM parameters are valid. */
  bool sanityCheck(const AtomTypesHolder& atomTypes) const;
  /** @brief Getter for the partial atomic charges for each atom*/
  std::vector<double> getChargesForEachAtom(const AtomTypesHolder& atomTypes) const;
  /**
   * @brief Resizes the C6 dispersion coefficient matrix to size NxN (N is number of distinct atom types).
   *        Furthermore, it creates a map which stores the indices of each atom type in the matrix.
   */
  void prepareC6Matrix(const AtomTypesHolder& atomTypes);
  /** @brief Getter for the C6 dispersion coefficient with indices of atom types given */
  float getC6(int indexOfAtomTypeA, int indexOfAtomTypeB) const;
  /** @brief Getter for the C6 dispersion coefficient for atoms with atom types a and b */
  float getC6(const std::string& a, const std::string& b) const;
  /** @brief Getter for the non-covalent parameters a1, s8, a2, beta and the atomic charges scaling factor */
  std::vector<double> getNonCovalentParameters() const;
  /** @brief Setter for the C6 dispersion coefficient with indices of atom types given */
  void setC6(int indexOfAtomTypeA, int indexOfAtomTypeB, float c6);
  /** @brief Setter for the C6 dispersion coefficient for atoms with atom types a and b */
  void setC6(const std::string& a, const std::string& b, float c6);
  /** @brief Resets the C6 matrix and the indices map */
  void resetC6Matrix();
  /** @brief Setter for the non-covalent parameters a1, s8, a2, beta and the atomic charges scaling factor */
  void setNonCovalentParameters(std::vector<double> nonCovalentParameters);
  /**
   * @brief Get MMDihedrals for the four atom types t1, t2, t3 and t4.
   *        (in principle there could be more than one (Fourier series))
   */
  std::vector<Dihedral> getMMDihedrals(std::string t1, std::string t2, std::string t3, std::string t4) const;
  /** @brief Get MMImproperDihedral for the four atom types t1, t2, t3 and t4 */
  ImproperDihedral getMMImproperDihedral(std::string central, std::string t2, std::string t3, std::string t4) const;
  /**
   * @brief Returns a const reference to the map containing the unique atom types of the system
   *        and their index in the c6 matrix.
   */
  const std::map<std::string, int>& getC6IndicesMap() const;

  // These functions add certain parameters of the MM model
  void addCharge(std::string atomType, double charge);
  void addDihedral(DihedralType dihedralType, DihedralParameters dihedralParameters);
  void addImproperDihedral(ImproperDihedralType improperDihedralType, ImproperDihedralParameters improperDihedralParameters);

  int evaluateNumDistinctAtomTypes(const AtomTypesHolder& atomTypes) const;

  // These functions return all parameters of a certain type by reference
  std::map<BondType, BondParameters>& getBonds();
  std::map<AngleType, AngleParameters>& getAngles();
  std::map<DihedralType, DihedralParameters>& getDihedrals();
  std::map<ImproperDihedralType, ImproperDihedralParameters>& getImproperDihedrals();
  std::map<std::string, double>& getCharges();

  // These functions return all parameters of a certain type by const. reference
  const std::map<BondType, BondParameters>& getBonds() const;
  const std::map<AngleType, AngleParameters>& getAngles() const;
  const std::map<DihedralType, DihedralParameters>& getDihedrals() const;
  const std::map<ImproperDihedralType, ImproperDihedralParameters>& getImproperDihedrals() const;
  const std::map<std::string, double>& getCharges() const;

 private:
  std::map<std::string, double> charges_;
  Eigen::MatrixXf c6Matrix_; // store as floats for efficiency
  std::map<std::string, int> c6IndicesMap_;
  std::vector<double> nonCovalentParameters_;
  std::map<DihedralType, DihedralParameters> dihedrals_;
  std::map<ImproperDihedralType, ImproperDihedralParameters> improperDihedrals_;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_SFAMPARAMETERS_H
