/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Swoose/Utilities/ConnectivityFileHandler.h>
#include <Swoose/Utilities/TopologyUtils.h>
#include <Utils/Bonds/BondOrderCollection.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace Scine;
using namespace SwooseUtilities;

void init_utilities(pybind11::module& m) {
  pybind11::class_<ConnectivityFileHandler> utilities(m, "utilities");

  utilities.def_static("write_connectivity_file", &ConnectivityFileHandler::writeListsOfNeighbors,
                       pybind11::arg("filename"), pybind11::arg("lists_of_neighbors"),
                       R"delim(
                         )delim");
}

void init_topology_utilities(pybind11::module& m) {
  pybind11::class_<TopologyUtils> topology_utilities(m, "topology_utilities");

  topology_utilities.def_static("generate_lists_of_neighbors",
                                &TopologyUtils::generateListsOfNeighborsFromBondOrderMatrix, pybind11::arg("natoms"),
                                pybind11::arg("bond_order_matrix"), pybind11::arg("minimal_bond_order_to_consider") = 0.5,
                                R"delim(
                         )delim");
}
