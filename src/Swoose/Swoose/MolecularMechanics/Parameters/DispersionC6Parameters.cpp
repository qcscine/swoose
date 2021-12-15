/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DispersionC6Parameters.h"
#include "../AtomTypesHolder.h"
#include "../SFAM/SfamParameters.h"
#include <Utils/Dftd3/Dftd3.h>
#include <Utils/Geometry/AtomCollection.h>

namespace Scine {
namespace MolecularMechanics {
namespace DispersionC6Parameters {

void fillC6MatrixForCurrentStructure(SfamParameters& parameters, const Utils::AtomCollection& structure,
                                     const AtomTypesHolder& atomTypes) {
  parameters.prepareC6Matrix(atomTypes);
  Utils::Dftd3::Dftd3 d3Method;
  std::vector<Utils::Dftd3::Dftd3Atom> structureOfDftd3Atoms;

  structureOfDftd3Atoms.reserve(structure.size());
  for (int index = 0; index < structure.size(); ++index) {
    Utils::Dftd3::Dftd3Atom newAtom(structure.getElement(index), structure.getPosition(index));
    newAtom.setIndex(index);
    structureOfDftd3Atoms.push_back(newAtom);
  }

  d3Method.setStructure(structureOfDftd3Atoms);

  for (auto& atom : structureOfDftd3Atoms) {
    double coordinationNumber = d3Method.calculateCoordinationNumber(atom);
    atom.setCoordinationNumber(coordinationNumber);
  }

  d3Method.setStructure(structureOfDftd3Atoms);

  const auto& indexMap = parameters.getC6IndicesMap();
  Eigen::MatrixXi occurrences(indexMap.size(), indexMap.size());
  occurrences.setZero();

  for (int a = 1; a < structure.size(); ++a) {
    for (int b = 0; b < a; ++b) {
      auto c6 = d3Method.calculateC6Coefficient(structureOfDftd3Atoms[a], structureOfDftd3Atoms[b]);
      auto atomTypeA = atomTypes.getAtomType(a);
      auto atomTypeB = atomTypes.getAtomType(b);
      int numOccurred = occurrences(indexMap.at(atomTypeA), indexMap.at(atomTypeB)) +
                        occurrences(indexMap.at(atomTypeB), indexMap.at(atomTypeA));
      if (numOccurred == 0) {
        parameters.setC6(atomTypeA, atomTypeB, static_cast<float>(c6));
      }
      else {
        float previousC6 = parameters.getC6(atomTypeA, atomTypeB);
        float newC6 = (previousC6 * numOccurred + c6) / (numOccurred + 1); // get the mean value
        parameters.setC6(atomTypeA, atomTypeB, newC6);
      }
      occurrences(indexMap.at(atomTypeA), indexMap.at(atomTypeB))++;
    }
  }
}

} // namespace DispersionC6Parameters
} // namespace MolecularMechanics
} // namespace Scine
