/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DihedralTerm.h"
#include "../MMExceptions.h"
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <Utils/Math/Tensor33.h>
#include <Eigen/Dense>

namespace Scine {
namespace MolecularMechanics {

DihedralTerm::DihedralTerm(AtomIndex firstAtom, AtomIndex secondAtom, AtomIndex thirdAtom, AtomIndex fourthAtom,
                           const Dihedral& dihedral, const DihedralType& typeOfDihedral)
  : firstAtom_(firstAtom),
    secondAtom_(secondAtom),
    thirdAtom_(thirdAtom),
    fourthAtom_(fourthAtom),
    dihedral_(dihedral),
    typeOfDihedral_(typeOfDihedral) {
}

DihedralTerm::~DihedralTerm() = default;

double DihedralTerm::evaluateDihedralTerm(const Utils::PositionCollection& positions,
                                          Utils::AtomicSecondDerivativeCollection& derivatives) const {
  if (this->disabled_)
    return 0.0;
  if (!dihedral_.hasParameters()) // Check this only if this term is not disabled
    throw MMDihedralParametersNotAvailableException(typeOfDihedral_.a1, typeOfDihedral_.a2, typeOfDihedral_.a3,
                                                    typeOfDihedral_.a4);

  // Calculating the dihedral energy and its derivatives according to
  // A. Blondel, M. Karplus, New Formulation for Derivatives of Torsion Angles and Improper Torsion Angles in Molecular
  // Mechanics: Elimination of Singularities, JCC, 17, 1996, 1132-1141.
  Eigen::Vector3d F(positions.row(firstAtom_) - positions.row(secondAtom_));
  Eigen::Vector3d G(positions.row(secondAtom_) - positions.row(thirdAtom_));
  Eigen::Vector3d H(positions.row(fourthAtom_) - positions.row(thirdAtom_));

  auto A = F.cross(G);
  auto B = H.cross(G);

  double A2 = A.squaredNorm();
  double B2 = B.squaredNorm();
  //  double A4 = A2 * A2;
  //  double B4 = B2 * B2;
  double G1 = G.norm();
  //  double G2 = G1 * G1;
  //  double G3 = G1 * G2;

  // Values needed for first derivatives
  auto gba2a = G1 / A2 * A;
  auto fgba2ga = F.dot(G) / (A2 * G1) * A;
  auto hgbb2gb = H.dot(G) / (B2 * G1) * B;
  auto gbb2b = G1 / B2 * B;

  // First derivatives:
  Eigen::Vector3d firstDer1 = -gba2a;
  Eigen::Vector3d firstDer2 = gba2a + fgba2ga - hgbb2gb;
  Eigen::Vector3d firstDer3 = hgbb2gb - fgba2ga - gbb2b;
  Eigen::Vector3d firstDer4 = gbb2b;

  //  // Values needed for second derivatives // TODO
  //  auto gca = G.cross(A);
  //  auto gcb = G.cross(B);
  //  auto fca = F.cross(A);
  //  auto hcb = H.cross(B);
  //  auto d2df2 = (Utils::tensor(A, gca) + Utils::tensor(gca, A)) * (G.norm() / A4);
  //  auto d2dh2 = (Utils::tensor(B, gcb) + Utils::tensor(gcb, B)) * (-G.norm() / B4);
  //  auto d2dg2 = (Utils::tensor(gca, A) + Utils::tensor(A, gca)) * (1.0 / (2 * G3 * A2)) +
  //               (Utils::tensor(A, fca) + Utils::tensor(fca, A)) * (F.dot(G) / (G.norm() * A4)) +
  //               (Utils::tensor(gcb, B) + Utils::tensor(B, gcb)) * (-1.0 / (2 * G3 * B2)) +
  //               (Utils::tensor(B, hcb) + Utils::tensor(hcb, B)) * (-H.dot(G) / (G.norm() * B4));
  //  auto d2dfdg = (Utils::tensor(-fca, A) * G.squaredNorm() + Utils::tensor(A, -gca) * F.dot(G)) * (1.0 / (G.norm() *
  //  A4)); auto d2dhdg = (Utils::tensor(-hcb, B) * G.squaredNorm() + Utils::tensor(B, -gcb) * H.dot(G)) * (-1.0 /
  //  (G.norm() * B4));

  // Derivative objects for theta
  Utils::AutomaticDifferentiation::Second3D h1, h2, h3, h4;

  //  // Add contributions to second derivatives // TODO
  //  setSecondDerivative(h1, d2df2, 1.0, 1.0);
  //  setSecondDerivative(h2, d2df2, -1.0, -1.0);
  //  setSecondDerivative(h2, d2dg2, 1.0, 1.0);
  //  setSecondDerivative(h2, d2dfdg, -1.0, 1.0);
  //  setSecondDerivative(h3, d2dh2, -1.0, -1.0);
  //  setSecondDerivative(h3, d2dg2, -1.0, -1.0);
  //  setSecondDerivative(h3, d2dhdg, -1.0, -1.0);
  //  setSecondDerivative(h4, d2dh2, 1.0, 1.0);

  // Add first derivatives contributions
  h1.setFirst3D(firstDer1);
  h2.setFirst3D(firstDer2);
  h3.setFirst3D(firstDer3);
  h4.setFirst3D(firstDer4);

  double theta = getTheta(A, B, G);
  auto result = dihedral_.getInteraction(theta);

  // Apply chain rule to get derivatives with respect to the energy
  derivatives[firstAtom_] += threeDimDer(result, h1);
  derivatives[secondAtom_] += threeDimDer(result, h2);
  derivatives[thirdAtom_] += threeDimDer(result, h3);
  derivatives[fourthAtom_] += threeDimDer(result, h4);

  return result.value();
}

Utils::AutomaticDifferentiation::Second3D DihedralTerm::threeDimDer(const Utils::AutomaticDifferentiation::Second1D& energy,
                                                                    const Utils::AutomaticDifferentiation::Second3D& alpha) const {
  return {energy.value(),
          energy.first() * alpha.dx(),
          energy.first() * alpha.dy(),
          energy.first() * alpha.dz(),
          energy.second() * alpha.dx() * alpha.dx() + energy.first() * alpha.XX(),
          energy.second() * alpha.dy() * alpha.dy() + energy.first() * alpha.YY(),
          energy.second() * alpha.dz() * alpha.dz() + energy.first() * alpha.ZZ(),
          energy.second() * alpha.dx() * alpha.dy() + energy.first() * alpha.XY(),
          energy.second() * alpha.dx() * alpha.dz() + energy.first() * alpha.XZ(),
          energy.second() * alpha.dy() * alpha.dz() + energy.first() * alpha.YZ()};
}

void DihedralTerm::setSecondDerivative(Utils::AutomaticDifferentiation::Second3D& h, const Utils::Tensor33& tensor,
                                       double derFirst, double derSecond) const {
  double f = derFirst * derSecond;
  h.setXX(h.XX() + f * tensor.x().x());
  h.setYY(h.YY() + f * tensor.y().y());
  h.setZZ(h.ZZ() + f * tensor.z().z());
  h.setXY(h.XY() + f * tensor.x().y() + f * tensor.y().x()); // TODO: DO I need 1/2 ?
  h.setXZ(h.XZ() + f * tensor.x().z() + f * tensor.z().x());
  h.setYZ(h.YZ() + f * tensor.y().z() + f * tensor.z().y());
}

double DihedralTerm::getTheta(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& G) {
  double acosArg = A.dot(B) / (A.norm() * B.norm());
  //  acosArg *= -1;
  double theta = acos(acosArg);
  // Needed because of numerical instabilities provoking theta = nan
  if (acosArg >= 1)
    theta = 0;
  else if (acosArg <= -1)
    theta = 4.0 * atan(1);

  // Invert sign of theta if it should be negative (NB: acos delivers only values between 0 and pi)
  double asinArg = B.cross(A).dot(G) / (A.norm() * B.norm() * G.norm());
  if (asinArg < 0) {
    theta *= -1;
  }

  return theta;
}

DihedralType DihedralTerm::getTypeOfDihedral() const {
  return typeOfDihedral_;
}

int DihedralTerm::getFirstAtom() const {
  return firstAtom_;
}

int DihedralTerm::getSecondAtom() const {
  return secondAtom_;
}

int DihedralTerm::getThirdAtom() const {
  return thirdAtom_;
}

int DihedralTerm::getFourthAtom() const {
  return fourthAtom_;
}

} // namespace MolecularMechanics
} // namespace Scine
