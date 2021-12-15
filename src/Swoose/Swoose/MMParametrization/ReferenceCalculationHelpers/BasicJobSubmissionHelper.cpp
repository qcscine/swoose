/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BasicJobSubmissionHelper.h"
#include "../MMParametrizationSettings.h"

namespace Scine {
namespace MMParametrization {
namespace BasicJobSubmissionHelper {

std::string determineMethodFamily(std::string method, std::string referenceProgram) {
  if (referenceProgram == SwooseUtilities::OptionNames::orcaOption ||
      referenceProgram == SwooseUtilities::OptionNames::turbomoleOption)
    return "dft";
  return method;
}

std::string determineBasisSet(const Utils::Settings& settings) {
  std::string referenceProgram = settings.getString(SwooseUtilities::SettingsNames::referenceProgram);
  if (referenceProgram == SwooseUtilities::OptionNames::sparrowOption || referenceProgram == SwooseUtilities::OptionNames::xtbOption)
    return "";
  return settings.getString(SwooseUtilities::SettingsNames::referenceBasisSet);
}

} // namespace BasicJobSubmissionHelper
} // namespace MMParametrization
} // namespace Scine
