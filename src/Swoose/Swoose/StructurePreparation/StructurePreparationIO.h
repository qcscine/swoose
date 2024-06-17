/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_STRUCTUREPREPARATIONIO_H
#define SWOOSE_STRUCTUREPREPARATIONIO_H

#include <boost/filesystem.hpp>
#include <string>
#include <vector>

namespace bfs = boost::filesystem;

namespace Scine {
namespace Core {
struct Log;
}
namespace StructurePreparation {
struct StructurePreparationData;
struct StructurePreparationFiles;
struct TitrableSite;

namespace StructurePreparationIO {
/**
 * @brief Converts an xyz file to an internal pdb file that fits the file format requirements for openbabel.
 */
void xyzToPdb(const std::string& xyzFile, const std::string& pdbFile);
/**
 * @brief Writes the atomic info to a file
 */
void writeAtomicInfoFileForProtein(const StructurePreparationData& data, const std::string& atomicInfoFile);
/**
 * @brief This function parses the user-generated atomic info file for the nonRegContainer, maps the nonRegContainer
 * indices to the full structure and appends the information to the global atomic info file.
 */
void addAtomicInformationForNonRegContainer(StructurePreparationFiles& files, std::vector<std::vector<int>> subsystemMapping);
// determines the suffix of a file
std::string getSuffix(const bfs::path& filepath);
/**
 * Writes the protein substructure to a PDB file.
 */
void writePdbFileWithResidueSpecifier(const StructurePreparationData& data, const std::string& proteinFile, Core::Log& log);

} // namespace StructurePreparationIO
} // namespace StructurePreparation
} // namespace Scine

#endif // SWOOSE_STRUCTUREPREPARATIONIO_H
