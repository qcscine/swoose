/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Swoose/MMParametrization/Parametrizer.h>
#include <Utils/Settings.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine;
using namespace MMParametrization;

void init_mm_parametrizer(pybind11::module& m) {
  pybind11::class_<Parametrizer> mm_parametrizer(m, "Parametrizer");

  mm_parametrizer.def(pybind11::init<>(),
                      R"delim(
                        Initialize the Parametrizer object.
                      )delim");

  mm_parametrizer.def(
      "parametrize_mm",
      [](Parametrizer& p, Utils::AtomCollection structure) -> void {
        pybind11::gil_scoped_release release;
        try {
          p.parametrize(structure);
        }
        catch (...) {
          pybind11::gil_scoped_acquire acquire;
          throw;
        }
        pybind11::gil_scoped_acquire acquire;
      },
      pybind11::arg("structure"),
      R"delim(
                           Perform an MM parametrization.
                           :param structure: The initial molecular structure for the simulation.
                         )delim");

  mm_parametrizer.def_property(
      "settings", [](Parametrizer& p) -> Scine::Utils::Settings& { return p.settings(); },
      [](Parametrizer& p, Scine::Utils::Settings settings) { p.settings() = std::move(settings); },
      pybind11::return_value_policy::reference, "Settings of the mm parametrizer.");

  mm_parametrizer.def_property(
      "log", [](Parametrizer& p) -> Scine::Core::Log& { return p.getLog(); },
      [](Parametrizer& p, Scine::Core::Log log) -> void { p.setLog(std::move(log)); },
      pybind11::return_value_policy::reference, "Logger of the MM parametrizer.");
}
