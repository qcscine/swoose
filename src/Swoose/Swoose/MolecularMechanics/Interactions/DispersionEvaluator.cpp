/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DispersionEvaluator.h"

namespace Scine {
namespace MolecularMechanics {

DispersionEvaluator::DispersionEvaluator(const Utils::AtomCollection& structure) : structure_(structure) {
}

double DispersionEvaluator::evaluate(Utils::FullSecondDerivativeCollection& derivatives,
                                     std::shared_ptr<Utils::Dftd3::Dftd3> d3, Eigen::MatrixXd& R0) {
  std::vector<Utils::Dftd3::Dftd3Atom> structureOfDftd3Atoms;

  structureOfDftd3Atoms.reserve(structure_.size());
  for (int index = 0; index < structure_.size(); ++index) {
    Utils::Dftd3::Dftd3Atom newAtom(structure_.getElement(index), structure_.getPosition(index));
    newAtom.setIndex(index);
    structureOfDftd3Atoms.push_back(newAtom);
  }

  d3->setStructure(structureOfDftd3Atoms);

  double energy = 0.0;

  for (auto& dispersion : dispersions_) {
    energy += dispersion.evaluateDispersionTerm(structureOfDftd3Atoms, derivatives, d3, R0);
  }

  return energy;
}

void DispersionEvaluator::setDispersionTerms(std::vector<DispersionTerm>&& dispersionTerms) {
  dispersions_ = dispersionTerms;
}

double DispersionEvaluator::a1_;
double DispersionEvaluator::s8_;
double DispersionEvaluator::a2_;
void DispersionEvaluator::setD3Parameters(std::vector<double> d3Parameters) {
  assert(d3Parameters.size() == 3);
  a1_ = d3Parameters.at(0);
  s8_ = d3Parameters.at(1);
  a2_ = d3Parameters.at(2);
}

double DispersionEvaluator::getA1() {
  return a1_;
}

double DispersionEvaluator::getS8() {
  return s8_;
}

double DispersionEvaluator::getA2() {
  return a2_;
}

} // namespace MolecularMechanics
} // namespace Scine
