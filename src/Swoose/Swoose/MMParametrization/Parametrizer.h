/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef MMPARAMETRIZATION_PARAMETRIZER_H
#define MMPARAMETRIZATION_PARAMETRIZER_H

#include "ParametrizationData.h"
#include <Core/Interfaces/MMParametrizer.h>

namespace Scine {

namespace Utils {
class Settings;
} // namespace Utils

namespace MMParametrization {
class ConnectivityGenerator;
class ReparametrizationHelper;

/**
 * @class Parametrizer Parametrizer.h
 * @brief This class manages the parametrization of an MM model from quantum-chemical reference data.
 */
class Parametrizer : public Core::MMParametrizer {
 public:
  static constexpr const char* model = "SFAM_parametrizer";
  /**
   * @brief Constructor.
   */
  Parametrizer();
  /**
   * @brief Main function of this class. It generates the MM parameters.
   * @param structure The entire system's structure.
   */
  void parametrize(Utils::AtomCollection structure) override;
  /**
   * @brief Accessor for the settings.
   */
  Utils::Settings& settings() override;
  /**
   * @brief Constant accessor for the settings.
   */
  const Utils::Settings& settings() const override;
  /**
   * @brief Getter for the name of the Parametrizer.
   */
  std::string name() const override;

 private:
  /*
   * @brief Performs some additional settings checks.
   */
  void performAdditionalSettingsChecks();
  /*
   * @brief Performs the initial set-up and fragmentation step.
   */
  void performInitialSetup(Utils::AtomCollection structure);
  /*
   * @brief Performs reference data generation step.
   */
  void generateReferenceData();
  /*
   * @brief Generates the topology and sets up the parameter optimization.
   */
  void setupParameterOptimization();
  /*
   * @brief Optimizes the parameters.
   */
  void optimizeParameters();
  /*
   * @brief Writes parameters and connectivity to files.
   */
  void writeParametersAndConnectivity();
  /*
   * @brief Generates the MM topology.
   */
  void generateTopology();
  /*
   * @brief Generates the atom types.
   */
  void generateAtomTypes();
  /*
   * @brief Sets the default settings values for method and basis set depending on the selected program.
   */
  void setDefaultsForMethodAndBasisSetSettings();
  /*
   * @brief Sets the protonation state of all pH sensitive sites after evaluation of the pKa.
   */
  void determineProtonationStateOfTitrableSites();
  // This object holds all the data used in the MM parametrization algorithm.
  ParametrizationData data_;
  // This object holds all the data used for titration.
  TitrationResults titrationResults_;
  // The settings.
  std::shared_ptr<Utils::Settings> settings_;
  // The connectivity generator class.
  std::shared_ptr<ConnectivityGenerator> connectivityGenerator_;
  // The re-parametrization helper class.
  std::shared_ptr<ReparametrizationHelper> reparametrizationHelper_;
};

} // namespace MMParametrization
} // namespace Scine

#endif // MMPARAMETRIZATION_PARAMETRIZER_H
