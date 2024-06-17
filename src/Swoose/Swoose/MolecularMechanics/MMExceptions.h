/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_MMEXCEPTIONS_H
#define MOLECULARMECHANICS_MMEXCEPTIONS_H

#include <exception>
#include <stdexcept>
#include <string>

/*
 * Defining own exceptions for molecular mechanics.
 */

namespace Scine {
namespace MolecularMechanics {

class MMBondParametersNotAvailableException : public std::exception {
 public:
  explicit MMBondParametersNotAvailableException(const std::string& atom1, const std::string& atom2) {
    message_ = "No bond parameters for " + atom1 + "-" + atom2;
  }

  const char* what() const noexcept override {
    return message_.c_str();
  }

 private:
  std::string message_;
};

class MMAngleParametersNotAvailableException : public std::exception {
 public:
  explicit MMAngleParametersNotAvailableException(const std::string& atom1, const std::string& atom2, const std::string& atom3) {
    message_ = "No angle parameters for " + atom1 + "-" + atom2 + "-" + atom3;
  }

  const char* what() const noexcept override {
    return message_.c_str();
  }

 private:
  std::string message_;
};

class MMDihedralParametersNotAvailableException : public std::exception {
 public:
  explicit MMDihedralParametersNotAvailableException(const std::string& atom1, const std::string& atom2,
                                                     const std::string& atom3, const std::string& atom4) {
    message_ = "No dihedral parameters for " + atom1 + "-" + atom2 + "-" + atom3 + "-" + atom4;
  }

  const char* what() const noexcept override {
    return message_.c_str();
  }

 private:
  std::string message_;
};

class MMImproperDihedralParametersNotAvailableException : public std::exception {
 public:
  explicit MMImproperDihedralParametersNotAvailableException(const std::string& atomCentral, const std::string& atom2,
                                                             const std::string& atom3, const std::string& atom4) {
    message_ = "No improper dihedral parameters for " + atomCentral + "-" + atom2 + "-" + atom3 + "-" + atom4;
  }

  const char* what() const noexcept override {
    return message_.c_str();
  }

 private:
  std::string message_;
};

class MMLjParametersNotAvailableException : public std::exception {
 public:
  explicit MMLjParametersNotAvailableException(const std::string& atomType) {
    message_ = "There are no Lennard Jones (vdW) parameters for " + atomType;
  }

  const char* what() const noexcept override {
    return message_.c_str();
  }

 private:
  std::string message_;
};

class ElementTypeIsNotAllowedForHydrogenBondsException : public std::runtime_error {
 public:
  explicit ElementTypeIsNotAllowedForHydrogenBondsException(const std::string& s) : runtime_error(s) {
  }
};

class TwoResultsForHydrogenBondsAreNotEqualException : public std::runtime_error {
 public:
  explicit TwoResultsForHydrogenBondsAreNotEqualException(const std::string& s) : runtime_error(s) {
  }
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_MMEXCEPTIONS_H
