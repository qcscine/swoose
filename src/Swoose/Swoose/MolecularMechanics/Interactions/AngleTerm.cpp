/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AngleTerm.h"
#include "../MMExceptions.h"
#include <Utils/Constants.h>
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Utils/Math/AutomaticDifferentiation/VectorDerivatives3D.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

namespace Scine {

using namespace Utils::AutomaticDifferentiation;
namespace MolecularMechanics {

AngleTerm::AngleTerm(AtomIndex firstAtom, AtomIndex secondAtom, AtomIndex thirdAtom, const Angle& angle, const AngleType& typeOfAngle)
  : firstAtom_(firstAtom), secondAtom_(secondAtom), thirdAtom_(thirdAtom), angle_(angle), typeOfAngle_(typeOfAngle) {
}

AngleTerm::~AngleTerm() = default;

double AngleTerm::evaluateAngleTerm(const Utils::PositionCollection& positions,
                                    Utils::AtomicSecondDerivativeCollection& derivatives) const {
  if (this->disabled_)
    return 0.0;
  if (!angle_.hasParameters()) // Check this only if this term is not disabled
    throw MMAngleParametersNotAvailableException(typeOfAngle_.a1, typeOfAngle_.a2, typeOfAngle_.a3);

  Eigen::Vector3d a(positions.row(firstAtom_) - positions.row(secondAtom_));
  Eigen::Vector3d b(positions.row(thirdAtom_) - positions.row(secondAtom_));

  double angle = acos((a.dot(b) / (a.norm() * b.norm())));
  Second1D result = angle_.getInteraction(angle);
  double energy = result.value();

  Second3D angleContributionAtom1;
  Second3D angleContributionAtom2;
  Second3D angleContributionAtom3;
  if (angle < Utils::Constants::pi - singularityCriterion_ && angle > singularityCriterion_) {
    calculateDerivativesWithNormalFormula(result, a, b, angleContributionAtom1, angleContributionAtom2, angleContributionAtom3);
  }
  else {
    calculateDerivativesForCriticalAngles(result, a, b, angleContributionAtom1, angleContributionAtom2, angleContributionAtom3);
    energy = 0;
  }

  derivatives[firstAtom_] += angleContributionAtom1;
  derivatives[secondAtom_] += angleContributionAtom2;
  derivatives[thirdAtom_] += angleContributionAtom3;

  return energy;
}

Second3D AngleTerm::threeDimDer(const Second1D& energy, const Second3D& alpha) {
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

void AngleTerm::calculateDerivativesWithNormalFormula(const Second1D& angleDerivative, const Eigen::Vector3d& a,
                                                      const Eigen::Vector3d& b, Second3D& angleContributionAtom1,
                                                      Second3D& angleContributionAtom2, Second3D& angleContributionAtom3) {
  // NB: suffix digit indicates with respect to which atom the derivative is taken
  auto a1 = VectorDerivatives3D::spatialVectorHessian3D(a);
  auto a2 = VectorDerivatives3D::spatialVectorHessian3DWithInverseDerivative(a);
  auto b3 = VectorDerivatives3D::spatialVectorHessian3D(b);
  auto b2 = VectorDerivatives3D::spatialVectorHessian3DWithInverseDerivative(b);

  Second3D alpha1 = arccos(a1.dot(b) / (a1.norm() * b.norm()));
  Second3D alpha2 = arccos(a2.dot(b2) / (a2.norm() * b2.norm())); // NB: maybe the second derivatives for atom 2 can be
                                                                  // obtained from atom 1 and 2 in the end and this is
                                                                  // not needed
  Second3D alpha3 = arccos(b3.dot(a) / (a.norm() * b3.norm()));

  angleContributionAtom1 = threeDimDer(angleDerivative, alpha1);
  angleContributionAtom2 = threeDimDer(angleDerivative, alpha2);
  angleContributionAtom3 = threeDimDer(angleDerivative, alpha3);
}

void AngleTerm::calculateDerivativesForCriticalAngles(const Second1D& angleDerivative, const Eigen::Vector3d& a,
                                                      const Eigen::Vector3d& b, Second3D& angleContributionAtom1,
                                                      Second3D& angleContributionAtom2, Second3D& angleContributionAtom3) {
  // Explanation of formulas: see Alain's derivation
  Eigen::Vector3d orientation = a.normalized();

  double inverseDistanceSum = 1.0 / a.norm() + 1.0 / b.norm();
  double d2thetad1 = 1.0 / a.squaredNorm();
  double d2thetad3 = 1.0 / b.squaredNorm();
  double d2thetad2 = inverseDistanceSum * inverseDistanceSum;

  double d2dxx1 = d2thetad1 * angleDerivative.second(); // Second derivative of the Energy with respect to x in local
                                                        // coordinate system for atom 1
  double d2dxx2 = d2thetad2 * angleDerivative.second();
  double d2dxx3 = d2thetad3 * angleDerivative.second();

  double x = orientation.x();
  double y = orientation.y();
  double z = orientation.z();

  // Get the local coordinate system: generate perpendicular unit x and y vectors
  Eigen::Vector3d yV(0, z, -y);
  if (std::abs(x) > std::abs(y))
    yV = Eigen::Vector3d(-z, 0, x);
  yV.normalize();
  Eigen::Vector3d xV = yV.cross(orientation);

  // Derivatives of local coordinates with respect to global coordinates
  double dxdx = xV.x();
  double dxdy = xV.y();
  double dxdz = xV.z();
  double dydx = yV.x();
  double dydy = yV.y();
  double dydz = yV.z();

  angleContributionAtom1 = Second3D(0, 0, 0, 0, d2dxx1 * (dxdx * dxdx + dydx * dydx), d2dxx1 * (dxdy * dxdy + dydy * dydy),
                                    d2dxx1 * (dxdz * dxdz + dydz * dydz), d2dxx1 * (dxdx * dxdy + dydx * dydy),
                                    d2dxx1 * (dxdx * dxdz + dydx * dydz), d2dxx1 * (dxdy * dxdz + dydy * dydz));
  angleContributionAtom2 = Second3D(0, 0, 0, 0, d2dxx2 * (dxdx * dxdx + dydx * dydx), d2dxx2 * (dxdy * dxdy + dydy * dydy),
                                    d2dxx2 * (dxdz * dxdz + dydz * dydz), d2dxx2 * (dxdx * dxdy + dydx * dydy),
                                    d2dxx2 * (dxdx * dxdz + dydx * dydz), d2dxx2 * (dxdy * dxdz + dydy * dydz));
  angleContributionAtom3 = Second3D(0, 0, 0, 0, d2dxx3 * (dxdx * dxdx + dydx * dydx), d2dxx3 * (dxdy * dxdy + dydy * dydy),
                                    d2dxx3 * (dxdz * dxdz + dydz * dydz), d2dxx3 * (dxdx * dxdy + dydx * dydy),
                                    d2dxx3 * (dxdx * dxdz + dydx * dydz), d2dxx3 * (dxdy * dxdz + dydy * dydz));
}

AngleType AngleTerm::getTypeOfAngle() const {
  return typeOfAngle_;
}

int AngleTerm::getFirstAtom() const {
  return firstAtom_;
}

int AngleTerm::getSecondAtom() const {
  return secondAtom_;
}

int AngleTerm::getThirdAtom() const {
  return thirdAtom_;
}

} // namespace MolecularMechanics
} // namespace Scine
