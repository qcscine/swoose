/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Files/tests_file_location.h"
#include <Core/Log.h>
#include <Swoose/MMParametrization/ParametrizationData.h>
#include <Swoose/Utilities/AtomicInformationReader.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
using namespace MMParametrization;
namespace Tests {

/**
 * @class AtomicInformationFileReaderTest AtomicInformationFileReaderTest.cpp
 * @brief Tests the AtomicInformationFileReader class
 * @test
 */
class AtomicInformationFileReaderTest : public Test {};

TEST_F(AtomicInformationFileReaderTest, AtomicInformationFileIsReadCorrectly) {
  ParametrizationData data;
  auto silentLogger = Core::Log::silent();
  SwooseUtilities::AtomicInformationReader reader(silentLogger);

  reader.read(atomic_information_file, data.formalCharges, data.unpairedElectrons, 450);

  ASSERT_THAT(data.formalCharges.size(), Eq(4));
  ASSERT_THAT(data.formalCharges[25], Eq(2));
  ASSERT_THAT(data.formalCharges[435], Eq(2));

  ASSERT_THAT(data.unpairedElectrons.size(), Eq(5));
  ASSERT_THAT(data.unpairedElectrons[25], Eq(1));
  ASSERT_THAT(data.unpairedElectrons[111], Eq(4));
}

TEST_F(AtomicInformationFileReaderTest, WrongAtomicInformationFileReadingFails) {
  ParametrizationData data;
  auto silentLogger = Core::Log::silent();
  SwooseUtilities::AtomicInformationReader reader(silentLogger);

  EXPECT_THROW(reader.read(atomic_information_file, data.formalCharges, data.unpairedElectrons, 400), std::runtime_error);
}

} // namespace Tests
} // namespace Scine