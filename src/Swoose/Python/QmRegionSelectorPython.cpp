/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Swoose/QMMM/QmRegionSelection/QmRegionSelector.h>
#include <Utils/Settings.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine;
using namespace Qmmm;

void init_qm_region_selector(pybind11::module& m) {
  pybind11::class_<QmRegionSelector> qm_region_selector(m, "QmRegionSelector");

  qm_region_selector.def(pybind11::init<>(),
                         R"delim(
                        Initialize the QmRegionSelector object.
                      )delim");
  qm_region_selector.def(
      "generate_qm_region",
      [](QmRegionSelector& s, const Utils::AtomCollection& structure) -> void {
        if (s.allowsPythonGILRelease()) {
          try {
            pybind11::gil_scoped_release release;
            s.generateQmRegion(structure);
          }
          catch (...) {
            pybind11::gil_scoped_acquire acquire;
            throw;
          }
          pybind11::gil_scoped_acquire acquire;
        }
        else {
          s.generateQmRegion(structure);
        }
      },
      pybind11::arg("structure"),
      R"delim(
                           Generates a (spherical) QM region around a central atom.
                           :param structure: The initial molecular structure.
                         )delim");
  qm_region_selector.def("set_underlying_calculator", &QmRegionSelector::setUnderlyingCalculator,
                         pybind11::arg("qmmmCalculator"),
                         R"delim(
                           Sets the underlying calculator for the QM region selection task.
                           :param structure: The QMMM calculator.
                         )delim");

  qm_region_selector.def("get_qm_region_structure", &QmRegionSelector::getQmRegionStructure,
                         R"delim(
                           Getter for the structure of the generated QM region.
                           :return The QM region structure.
                         )delim");

  qm_region_selector.def("get_qm_region_indices", &QmRegionSelector::getQmRegionIndicesWithoutLinkAtoms,
                         R"delim(
                           Getter for the atom indices of the generated QM region.
                           :return A list of QM region indices.
                         )delim");

  qm_region_selector.def("get_qm_region_charge_multiplicity", &QmRegionSelector::getQmRegionChargeAndMultiplicity,
                         R"delim(
                           Getter for the atom indices of the generated QM region.
                           :return A pair of charge and multiplicity for the QM region.
                         )delim");

  qm_region_selector.def_property(
      "settings", [](QmRegionSelector& s) -> Scine::Utils::Settings& { return s.settings(); },
      [](QmRegionSelector& s, Scine::Utils::Settings settings) { s.settings() = std::move(settings); },
      pybind11::return_value_policy::reference, "Settings of the qm region selector.");
}
