/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ParameterFileWriter.h"
#include <Swoose/MolecularMechanics/SFAM/SfamParameters.h>
#include <fstream>

namespace Scine {
namespace MMParametrization {

void ParameterFileWriter::writeSfamParametersToFile(const std::string& filename,
                                                    const MolecularMechanics::SfamParameters& parameters) {
  std::ofstream parFile(filename);

  parFile << "# MM parameters of SFAM generated by SCINE\n\n";

  parFile << "! bonds" << std::endl;
  if (parameters.getBonds().empty())
    parFile << "*" << std::endl;
  for (const auto& bond : parameters.getBonds()) {
    parFile << bond.first.a1 << "   " << bond.first.a2 << "   " << bond.second.getEquilibriumBondLength() << "   "
            << bond.second.getForceConstant() << std::endl;
  }

  parFile << "\n";

  parFile << "! angles" << std::endl;
  if (parameters.getAngles().empty())
    parFile << "*" << std::endl;
  for (const auto& angle : parameters.getAngles()) {
    parFile << angle.first.a1 << "   " << angle.first.a2 << "   " << angle.first.a3 << "   "
            << angle.second.getEquilibriumAngle() << "   " << angle.second.getForceConstant() << std::endl;
  }

  parFile << "\n";

  parFile << "! dihedrals" << std::endl;
  if (parameters.getDihedrals().empty())
    parFile << "*" << std::endl;
  for (const auto& dihedral : parameters.getDihedrals()) {
    parFile << dihedral.first.a1 << "   " << dihedral.first.a2 << "   " << dihedral.first.a3 << "   "
            << dihedral.first.a4 << "   " << dihedral.second.getHalfBarrierHeight() << "   "
            << dihedral.second.getPhaseShift() << "   " << dihedral.second.getPeriodicity() << std::endl;
  }

  parFile << "\n";

  parFile << "! impropers" << std::endl;
  if (parameters.getImproperDihedrals().empty())
    parFile << "*" << std::endl;
  for (const auto& improperDihedral : parameters.getImproperDihedrals()) {
    parFile << improperDihedral.first.ac << "   " << improperDihedral.first.a2 << "   " << improperDihedral.first.a3
            << "   " << improperDihedral.first.a4 << "   " << improperDihedral.second.getEquilibriumAngle() << "   "
            << improperDihedral.second.getForceConstant() << std::endl;
  }

  parFile << "\n";

  parFile << "! charges" << std::endl;
  if (parameters.getCharges().empty())
    parFile << "*" << std::endl;
  for (const auto& charge : parameters.getCharges()) {
    parFile << charge.first << "   " << charge.second << std::endl;
  }

  parFile << "\n";

  parFile << "! non-covalent" << std::endl;
  auto params = parameters.getNonCovalentParameters();
  if (params.empty())
    parFile << "*" << std::endl;
  parFile << "a1"
          << "   " << params.at(0) << std::endl;
  parFile << "s8"
          << "   " << params.at(1) << std::endl;
  parFile << "a2"
          << "   " << params.at(2) << std::endl;
  parFile << "beta"
          << "   " << params.at(3) << std::endl;
  parFile << "scaling_factor_for_atomic_charges"
          << "   " << params.at(4) << std::endl;

  parFile << "\n";

  parFile << "! c6 coefficients" << std::endl;
  for (int a = 0; a < parameters.getC6IndicesMap().size(); ++a) {
    for (int b = 0; b <= a; ++b) {
      auto c6 = parameters.getC6(a, b);
      parFile << c6 << "  ";
    }
    parFile << "\n";
  }
}

} // namespace MMParametrization
} // namespace Scine