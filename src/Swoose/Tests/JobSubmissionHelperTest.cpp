/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Swoose/MMParametrization/MMParametrizationSettings.h>
#include <Swoose/MMParametrization/Parametrizer.h>
#include <gmock/gmock.h>
#include <Swoose/MMParametrization/ReferenceCalculationHelpers/BasicJobSubmissionHelper.cpp>

using namespace testing;
namespace Scine {
using namespace MMParametrization;
using namespace BasicJobSubmissionHelper;
namespace Tests {

/**
 * @class JobSubmissionHelperTest JobSubmissionHelperTest.cpp
 * @brief Tests helper function for job submission in MM parametrization.
 * @test
 */
class JobSubmissionHelperTest : public Test {};

TEST_F(JobSubmissionHelperTest, MethodFamilyIsCorrectlyDeterminedFromMethodAndReferenceProgram) {
  auto result = determineMethodFamily("pbe d3bj", "orca");
  ASSERT_STREQ(result.c_str(), "dft");
  result = determineMethodFamily("pm3", "orca");
  ASSERT_STREQ(result.c_str(), "dft");
  result = determineMethodFamily("pbe", "turbomole");
  ASSERT_STREQ(result.c_str(), "dft");
  result = determineMethodFamily("gfn2", "xtb");
  ASSERT_STREQ(result.c_str(), "gfn2");
  result = determineMethodFamily("pm6", "sparrow");
  ASSERT_STREQ(result.c_str(), "pm6");
  result = determineMethodFamily("dftb3", "sparrow");
  ASSERT_STREQ(result.c_str(), "dftb3");
}

TEST_F(JobSubmissionHelperTest, BasisSetIsCorrectlyDeterminedFromSettings) {
  Parametrizer parametrizer;
  std::string arbitraryBasisSetString = "any basis";
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceBasisSet, arbitraryBasisSetString);

  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceProgram, "orca");
  ASSERT_STREQ(determineBasisSet(parametrizer.settings()).c_str(), arbitraryBasisSetString.c_str());
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceProgram, "turbomole");
  ASSERT_STREQ(determineBasisSet(parametrizer.settings()).c_str(), arbitraryBasisSetString.c_str());
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceProgram, "xtb");
  ASSERT_STREQ(determineBasisSet(parametrizer.settings()).c_str(), "");
  parametrizer.settings().modifyString(SwooseUtilities::SettingsNames::referenceProgram, "sparrow");
  ASSERT_STREQ(determineBasisSet(parametrizer.settings()).c_str(), "");
}

} // namespace Tests
} // namespace Scine
