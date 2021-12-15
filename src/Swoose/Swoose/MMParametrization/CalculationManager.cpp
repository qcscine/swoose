/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CalculationManager.h"
#include "MMParametrizationSettings.h"
#include "ParametrizationData.h"
#include "ReferenceCalculationHelpers/DirectCalculationsHelper.h"
#include "ReferenceCalculationHelpers/ReferenceCalculationsIO.h"
#include <Core/Log.h>
#include <numeric>

#ifdef SWOOSE_COMPILE_DATABASE
#  include "ReferenceCalculationHelpers/DatabaseHelper.h"
#endif

namespace Scine {
namespace MMParametrization {

CalculationManager::CalculationManager(ParametrizationData& data, std::shared_ptr<Utils::Settings> settings, Core::Log& log)
  : data_(data), settings_(settings), log_(log) {
  baseWorkingDirectory_ = settings_->getString(Utils::ExternalQC::SettingsNames::baseWorkingDirectory);
  referenceDataDirectory_ = settings->getString(SwooseUtilities::SettingsNames::referenceDataDirectory);

  // Get mode from settings
  mode_ = settings_->getString(SwooseUtilities::SettingsNames::referenceDataMode);

  // Get reference program from settings
  referenceProgram_ = settings_->getString(SwooseUtilities::SettingsNames::referenceProgram);
}

void CalculationManager::calculateReferenceData() {
  using namespace SwooseUtilities::OptionNames;
  auto numberOfStructures = data_.vectorOfStructures.size();

  // Resize the vector of Hessian matrices
  data_.vectorOfHessians.resize(numberOfStructures);

  if (mode_ == writeToFilesMode) {
    log_.output << "Writing data for all molecular structures to directory: " << referenceDataDirectory_ << Core::Log::endl;
    ReferenceCalculationsIO::writeXyzFiles(data_, referenceDataDirectory_);
    log_.output << "Done." << Core::Log::endl << Core::Log::endl;
  }
  else if (mode_ == readFromFilesMode) {
    log_.output << "Reading reference data from disk..." << Core::Log::endl;
    ReferenceCalculationsIO::readReferenceDataFromFiles(data_, referenceDataDirectory_, settings_, log_);
    log_.output << "Done." << Core::Log::endl << Core::Log::endl;
  }
  else if (mode_ == directMode) {
    DirectCalculationsHelper::performReferenceCalculations(data_, settings_, baseWorkingDirectory_, log_);
    if (settings_->getBool(SwooseUtilities::SettingsNames::useGaussianOptionKey))
      DirectCalculationsHelper::calculateAtomicChargesWithGaussian(data_, settings_, baseWorkingDirectory_, log_);
  }
  else if (mode_ == databaseMode) {
#ifdef SWOOSE_COMPILE_DATABASE
    DatabaseHelper db(data_, settings_, log_);
    db.runCalculationsAndCollectResults();
#else
    throw std::runtime_error("Swoose was not compiled with database support.");
#endif
  }
  double sumOfAtomicCharges = std::accumulate(data_.atomicCharges.begin(), data_.atomicCharges.end(), 0.0);
  log_.debug << "Sum of all atomic partial charges: " << sumOfAtomicCharges << Core::Log::endl;
}

} // namespace MMParametrization
} // namespace Scine
