/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BasicJobSubmissionHelper.h"
#include "OptionNames.h"
#include "SettingsNames.h"

namespace Scine {
namespace SwooseUtilities {
namespace BasicJobSubmissionHelper {

std::string determineMethodFamily(std::string method, std::string referenceProgram) {
  if (referenceProgram == OptionNames::orcaOption || referenceProgram == OptionNames::turbomoleOption)
    return "dft";
  return method;
}

std::string determineBasisSet(const Utils::Settings& settings) {
  std::string referenceProgram = settings.getString(SettingsNames::referenceProgram);
  if (referenceProgram == OptionNames::sparrowOption || referenceProgram == OptionNames::xtbOption)
    return "";
  return settings.getString(SettingsNames::referenceBasisSet);
}

} // namespace BasicJobSubmissionHelper
} // namespace SwooseUtilities
} // namespace Scine
