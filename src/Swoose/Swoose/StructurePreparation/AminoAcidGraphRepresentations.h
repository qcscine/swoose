/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef PDBPREPARATION_AMINOACIDGRAPHREPRESENTATIONS_H
#define PDBPREPARATION_AMINOACIDGRAPHREPRESENTATIONS_H

#include <array>
#include <map>
#include <string>
#include <vector>

namespace Scine {

namespace Molassembler {
class Graph;
}

namespace StructurePreparation {
namespace AminoAcids {
/**
 * @brief Gets a map of amino acid names and corresponding molassembler graphs.
 * @return std::map<std::string, Molassembler::Graph>
 */
std::map<std::string, Molassembler::Graph> getAminoAcidGraphs();
/**
 * @brief Returns a vector of atom types for a specific amino acid.
 */
std::vector<std::string> getResidueTypes(const std::string& residueName);
// A vector of amino acids which should be queried in the nanoscale structure in this order.
static constexpr std::array<const char*, 24> aminoAcidHierarchy = {
    "PYL", "PRO", "HIS", "TRP", "TYR", "PHE", "LYS", "LEU", "ILE", "VAL", "GLU", "GLN",
    "ARG", "ASN", "ASP", "MSE", "MET", "THR", "SER", "DCY", "SEC", "CYS", "ALA", "GLY"};

} // namespace AminoAcids
} // namespace StructurePreparation
} // namespace Scine

#endif // PDBPREPARATION_AMINOACIDGRAPHREPRESENTATIONS_H
