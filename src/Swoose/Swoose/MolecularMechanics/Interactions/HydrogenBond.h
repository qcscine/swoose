/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MOLECULARMECHANICS_HYDROGENBOND_H
#define MOLECULARMECHANICS_HYDROGENBOND_H

#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Math/AutomaticDifferentiation/Second1D.h>
#include <array>

namespace Scine {
namespace MolecularMechanics {
namespace HydrogenBondHelper {

struct Donor {
  Donor(std::string hAtomType, std::string dAtomType, double kHBD, double chargeD);
  bool operator==(const Donor& rhs) const;
  std::string hAtomType, dAtomType;
  double kHBD, chargeD;
};

inline Donor::Donor(std::string dAtomType, std::string hAtomType, double kHBD, double chargeD)
  : hAtomType(std::move(hAtomType)), dAtomType(std::move(dAtomType)), kHBD(std::move(kHBD)), chargeD(std::move(chargeD)) {
}

inline bool Donor::operator==(const Donor& rhs) const {
  if ((dAtomType == rhs.dAtomType) && (hAtomType == rhs.hAtomType) && (kHBD == rhs.kHBD) && (chargeD == rhs.chargeD)) {
    return true;
  }
  else {
    return false;
  }
}

struct Acceptor {
  Acceptor(std::string aAtomType, double kHBA, double chargeA);
  bool operator==(const Acceptor& rhs) const;
  std::string aAtomType;
  double kHBA, chargeA;
};

inline Acceptor::Acceptor(std::string aAtomType, double kHBA, double chargeA)
  : aAtomType(std::move(aAtomType)), kHBA(std::move(kHBA)), chargeA(std::move(chargeA)) {
}

inline bool Acceptor::operator==(const Acceptor& rhs) const {
  if ((aAtomType == rhs.aAtomType) && (kHBA == rhs.kHBA) && (chargeA == rhs.chargeA)) {
    return true;
  }
  else {
    return false;
  }
}

static constexpr std::array<Utils::ElementType, 4> vectorOfDonorOrAcceptorElements_ = {
    Utils::ElementType::N, Utils::ElementType::O, Utils::ElementType::F, Utils::ElementType::Cl};
} // namespace HydrogenBondHelper

/**
 * @class HydrogenBond HydrogenBond.h
 * @brief Class calculating the energy and derivatives
 *        for a hydrogen bond based solely on the distance or on the angle, i.e. in 1 dimension, respectively.
 */
class HydrogenBond {
 public:
  /** @brief Default Constructor. */
  HydrogenBond();

  /**
   * @brief Calculate energy contribution from the distance with derivatives.
   */
  Utils::AutomaticDifferentiation::Second1D getInteractionDistanceVariable(double distance, double angle, double chargeDonor,
                                                                           double chargeAcceptor, double constantDonor,
                                                                           double constantAcceptor) const;
  /**
   * @brief Calculate energy contribution from
   * the angle with derivatives.
   */
  Utils::AutomaticDifferentiation::Second1D getInteractionAngleVariable(double distance, double angle, double chargeDonor,
                                                                        double chargeAcceptor, double constantDonor,
                                                                        double constantAcceptor) const;
};

} // namespace MolecularMechanics
} // namespace Scine

#endif // MOLECULARMECHANICS_HYDROGENBOND_H
