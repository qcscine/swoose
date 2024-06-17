/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <pybind11/pybind11.h>

void init_tasks(pybind11::module& m);
void init_mm_parametrizer(pybind11::module& m);
void init_qm_region_selector(pybind11::module& m);
void init_utilities(pybind11::module& m);
void init_topology_utilities(pybind11::module& m);

PYBIND11_MODULE(scine_swoose, m) {
  m.doc() = "Pybind11 Bindings for SCINE Swoose";
  init_tasks(m);
  init_mm_parametrizer(m);
  init_qm_region_selector(m);
  init_utilities(m);
  init_topology_utilities(m);
}