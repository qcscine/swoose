/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GaffParameterDefaultsProvider.h"
#include "../AtomTypesHolder.h"
#include "../Parameters/AngleParameters.h"
#include "../Parameters/BondParameters.h"
#include "../Parameters/LennardJonesParameters.h"
#include "GaffDefaultParameters.h"
#include "GaffParameters.h"

namespace Scine {
namespace MolecularMechanics {

std::unique_ptr<GaffParameters> GaffParameterDefaultsProvider::getParameters() {
  auto parameters = std::make_unique<GaffParameters>();
  processParameters(*parameters);
  return parameters;
}

void GaffParameterDefaultsProvider::processParameters(GaffParameters& parameters) {
  using namespace GaffDefaultParameters;

  assert(bondAtomTypes.size() == bondParameters.size());
  for (int i = 0; i < bondAtomTypes.size(); ++i) {
    // In GAFF, the harmonic potential is k*x^2 instead of 0.5*k*x^2. Hence, the "2.0 * ".
    parameters.addBond(BondType(bondAtomTypes[i][0], bondAtomTypes[i][1]),
                       BondParameters(2.0 * bondParameters[i][0], bondParameters[i][1]));
  }

  assert(anglesAtomTypes.size() == anglesParameters.size());
  for (int i = 0; i < anglesAtomTypes.size(); ++i) {
    // In GAFF, the harmonic potential is k*x^2 instead of 0.5*k*x^2. Hence, the "2.0 * ".
    parameters.addAngle(AngleType(anglesAtomTypes[i][0], anglesAtomTypes[i][1], anglesAtomTypes[i][2]),
                        AngleParameters(2.0 * anglesParameters[i][0], anglesParameters[i][1]));
  }

  assert(dihedralsAtomTypes.size() == dihedralsParameters.size());
  for (int i = 0; i < dihedralsAtomTypes.size(); ++i) {
    const double halfBarrier = dihedralsParameters[i][1] / dihedralsParameters[i][0];
    parameters.addDihedral(DihedralType(dihedralsAtomTypes[i][0], dihedralsAtomTypes[i][1], dihedralsAtomTypes[i][2],
                                        dihedralsAtomTypes[i][3]),
                           DihedralParameters(halfBarrier, dihedralsParameters[i][2], dihedralsParameters[i][3]));
  }

  assert(improperDihedralsAtomTypes.size() == improperDihedralsParameters.size());
  for (int i = 0; i < improperDihedralsAtomTypes.size(); ++i) {
    parameters.addImproperDihedral(ImproperDihedralType(improperDihedralsAtomTypes[i][2], improperDihedralsAtomTypes[i][0],
                                                        improperDihedralsAtomTypes[i][1], improperDihedralsAtomTypes[i][3]),
                                   DihedralParameters(improperDihedralsParameters[i][0], improperDihedralsParameters[i][1],
                                                      improperDihedralsParameters[i][2]));
  }

  assert(ljAtomTypes.size() == ljParameters.size());
  for (int i = 0; i < ljAtomTypes.size(); ++i) {
    parameters.addLennardJones(ljAtomTypes[i], LennardJonesParameters(ljParameters[i][0], ljParameters[i][1]));
  }
}

} // namespace MolecularMechanics
} // namespace Scine
