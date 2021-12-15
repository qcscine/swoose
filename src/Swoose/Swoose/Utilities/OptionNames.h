/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_OPTIONNAMES_H
#define SWOOSEUTILITIES_OPTIONNAMES_H

#include <Utils/Settings.h>

namespace Scine {
namespace SwooseUtilities {
namespace OptionNames {

static constexpr const char* redistributedChargeOption = "rc";
static constexpr const char* redistributedChargeAndDipolesOption = "rcd";
static constexpr const char* directMode = "direct";
static constexpr const char* databaseMode = "database";
static constexpr const char* writeToFilesMode = "write";
static constexpr const char* readFromFilesMode = "read";
static constexpr const char* orcaOption = "orca";
static constexpr const char* turbomoleOption = "turbomole";
static constexpr const char* sparrowOption = "sparrow";
static constexpr const char* xtbOption = "xtb";

} // namespace OptionNames
} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_OPTIONNAMES_H
