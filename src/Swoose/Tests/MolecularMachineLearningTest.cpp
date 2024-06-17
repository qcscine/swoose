/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TestUtilities/ReferenceDataForTests.h"
#include <Swoose/MachineLearning/MolecularMachineLearningModel.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/MolecularTrajectory.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
using namespace Swoose::MachineLearning;
namespace Tests {

/**
 * @class AMolecularMachineLearningTest MolecularMachineLearningTest.cpp
 * @brief Tests the molecular machine learning class to predict energies and atomic forces from reference data.
 * @test
 */
class AMolecularMachineLearningTest : public Test {
 public:
  Utils::MolecularTrajectory trajectory;
  Utils::MolecularTrajectory smallTrajectory;
  std::vector<double> energies;
  std::vector<Swoose::ForcesCollection> refForces;

  void SetUp() override {
    trajectory = Utils::MolecularTrajectoryIO::read(Utils::MolecularTrajectoryIO::format::xyz, butane_trajectory_xyz_file);

    smallTrajectory.setElementTypes(trajectory.getElementTypes());
    for (int i = 0; i < 100; ++i) {
      smallTrajectory.push_back(trajectory[i]);
    }
  }
};

TEST_F(AMolecularMachineLearningTest, ModelLearnsEnergiesAccurately) {
  bool debug = true;

#ifdef NDEBUG
  debug = false;
#endif

  std::vector<double> energies = Swoose::ReferenceDataForTests::parseReferenceEnergies();
  std::vector<double> allEnergies = energies;
  energies.resize(100);

  smallTrajectory.setEnergies(energies);
  trajectory.setEnergies(allEnergies);
  std::vector<Swoose::ForcesCollection> emptyForces;

  MolecularMachineLearningModel model;
  Eigen::VectorXd sigma(1);
  sigma[0] = 50.0;
  model.energyPredictor().setKernel(Utils::MachineLearning::Kernels::laplacianKernel, sigma);

  model.setReferenceData(smallTrajectory, emptyForces);

  auto result = model.evaluateEnergyModel(5);

  ASSERT_TRUE(result.first < 2.6);
  ASSERT_TRUE(result.first > 2.2);
  ASSERT_TRUE(result.second < 0.8);
  ASSERT_TRUE(result.second > 0.4);

  model.setReferenceData(trajectory, emptyForces);

  if (!debug) {
    result = model.evaluateEnergyModel(5);
    ASSERT_TRUE(result.first < 1.4);
    ASSERT_TRUE(result.second < 0.1);
  }

  std::stringstream ss("14\n\n"
                       "C     -0.1517065048   -0.6961577805    1.8429090617\n"
                       "H     -1.0605159877   -0.0912403196    2.1350466010\n"
                       "C     -0.0322314107   -0.7154656646    0.2985738962\n"
                       "C     -1.2932067816   -1.4688657912   -0.3274061427\n"
                       "C     -1.3698090055   -1.4046089176   -1.8253979208\n"
                       "H     -0.3789796045   -1.6782849470    2.1597735031\n"
                       "H      0.7475565292   -0.3719165650    2.4185762054\n"
                       "H      0.9583765228   -1.0825957352   -0.0678636930\n"
                       "H      0.0265105106    0.3245839593   -0.0176375150\n"
                       "H     -2.2462868113   -1.0580789863    0.2042904568\n"
                       "H     -1.3967805910   -2.5245723804    0.0491738773\n"
                       "H     -2.0826368896   -2.1956554317   -2.1709125491\n"
                       "H     -0.4564330954   -1.6136982029   -2.4058824289\n"
                       "H     -1.7212601432   -0.4012767498   -2.1696525501\n");

  auto newStructure = Utils::XyzStreamHandler::read(ss);
  model.trainEnergyModel();
  auto predictedEnergy = model.predictEnergy(newStructure);
  ASSERT_THAT(predictedEnergy, DoubleNear(10.3784, 0.5)); // This structure's energy should be predicted rather well.

  // With different kernel (5th degree polynomial)
  Eigen::VectorXd hyperparams(3);
  hyperparams << 5.0, 0.1, 1.0;
  model.energyPredictor().setKernel(Utils::MachineLearning::Kernels::polynomialKernel, hyperparams);
  model.trainEnergyModel();
  predictedEnergy = model.predictEnergy(newStructure);
  ASSERT_THAT(predictedEnergy, DoubleNear(10.3784, 1.0)); // Does not predict this structure as accurately.

  // With linear kernel
  model.energyPredictor().setKernel(Utils::MachineLearning::Kernels::linearKernel, {});
  model.trainEnergyModel();
  predictedEnergy = model.predictEnergy(newStructure);
  ASSERT_THAT(predictedEnergy, DoubleNear(10.3784, 2.5)); // Rather inaccurate, but still under 2.5 kcal/mol deviation.
}

TEST_F(AMolecularMachineLearningTest, ModelLearnsForcesAccurately) {
  bool debug = true;

#ifdef NDEBUG
  debug = false;
#endif

  refForces = Swoose::ReferenceDataForTests::parseReferenceForces(reference_forces_butane_file, 1000, 14);

  auto allRefForces = refForces;
  refForces.resize(100);

  std::stringstream ss("14\n\n"
                       "C     -0.1517065048   -0.6961577805    1.8429090617\n"
                       "H     -1.0605159877   -0.0912403196    2.1350466010\n"
                       "C     -0.0322314107   -0.7154656646    0.2985738962\n"
                       "C     -1.2932067816   -1.4688657912   -0.3274061427\n"
                       "C     -1.3698090055   -1.4046089176   -1.8253979208\n"
                       "H     -0.3789796045   -1.6782849470    2.1597735031\n"
                       "H      0.7475565292   -0.3719165650    2.4185762054\n"
                       "H      0.9583765228   -1.0825957352   -0.0678636930\n"
                       "H      0.0265105106    0.3245839593   -0.0176375150\n"
                       "H     -2.2462868113   -1.0580789863    0.2042904568\n"
                       "H     -1.3967805910   -2.5245723804    0.0491738773\n"
                       "H     -2.0826368896   -2.1956554317   -2.1709125491\n"
                       "H     -0.4564330954   -1.6136982029   -2.4058824289\n"
                       "H     -1.7212601432   -0.4012767498   -2.1696525501\n");

  auto newStructure = Utils::XyzStreamHandler::read(ss);

  MolecularMachineLearningModel model;

  if (debug)
    model.setReferenceData(smallTrajectory, refForces);
  else
    model.setReferenceData(trajectory, allRefForces);

  Eigen::VectorXd hyperparams(3);
  hyperparams << 5.0, 0.1, 1.0;
  for (int i = 0; i < trajectory.molecularSize(); ++i) {
    model.forcePredictor(i).setKernel(Utils::MachineLearning::Kernels::polynomialKernel, hyperparams);
  }

  // 5-fold cross validation with pooled variance
  auto result = model.evaluateForcesModel(5);

  if (debug) {
    ASSERT_TRUE(result.first < 1.6);
    ASSERT_TRUE(result.first > 1.3);
    ASSERT_TRUE(result.second > 0.1);
    ASSERT_TRUE(result.second < 0.3);
  }
  else {
    ASSERT_TRUE(result.first < 1.0);
    ASSERT_TRUE(result.second < 0.05);
  }

  auto oldMeanAbsoluteError = result.first;

  // 5-fold cross validation with true combined variance
  result = model.evaluateForcesModel(5, false);
  ASSERT_THAT(oldMeanAbsoluteError, DoubleNear(result.first, 1e-4));

  if (debug) {
    ASSERT_TRUE(result.second > 1.0);
    ASSERT_TRUE(result.second < 1.3);
  }
  else {
    ASSERT_TRUE(result.second > 0.05);
    ASSERT_TRUE(result.second < 1.0);
  }

  model.trainForcesModel();
  auto predictedForces = model.predictForces(newStructure);
  ASSERT_THAT(predictedForces.rows(), Eq(14));
  ASSERT_THAT(predictedForces.cols(), Eq(3));

  if (!debug) {
    Eigen::MatrixXd trueGradient(14, 3);
    trueGradient << 0.74334696, -8.6716958, 4.7221861, -4.3623572, 0.66275233, 0.37364911, 4.6895998, 4.2817808,
        -3.7440435, -1.887762, -4.4013384, -12.39901, 7.1070462, 0.6737955, 4.9976042, 1.8836753, 12.10442, -4.9575571,
        1.3633757, -1.1950653, 2.1242209, 2.2497201, 3.7596078, -1.1497938, 2.8652909, -4.924923, 3.471784, -6.9714376,
        -0.15802766, 6.3884607, -4.9618112, -1.2622967, 1.4547107, 0.43898267, -3.0377511, 0.87246102, -3.3617611,
        0.35965852, -2.2593211, 0.20409146, 1.8090834, 0.10464879;

    Eigen::MatrixXd trueForces = -trueGradient; // convert to forces

    int accuratePredictions = 0;
    for (int n = 0; n < 14; ++n) {
      for (int m = 0; m < 3; ++m) {
        if (std::abs(trueForces(n, m) - predictedForces(n, m)) < 1.0)
          accuratePredictions++;
      }
    }

    ASSERT_TRUE(accuratePredictions > 34); // more than 80 percent should be accurate predictions
  }
}

} // namespace Tests
} // namespace Scine
