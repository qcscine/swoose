/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_QMMM_QMREGIONCANDIDATEGENERATOR_H
#define SWOOSE_QMMM_QMREGIONCANDIDATEGENERATOR_H

#include <list>
#include <vector>

namespace Scine {

namespace Core {
class Log;
} // namespace Core

namespace Utils {
class AtomCollection;
class BondOrderCollection;
class Settings;
} // namespace Utils

namespace SwooseUtilities {
class SubsystemGenerator;
class FragmentAnalyzer;
} // namespace SwooseUtilities

namespace Qmmm {
struct QmmmModel;
namespace QmRegionCandidateGenerator {

/**
 * @brief Generates QM Region candidates with given settings and stores these in the vectors qmRegionCandidates
 *        and qmRegionCandidateIndices.
 * @param qmmmModelCandidates A vector that stores the QM/MM candidate models.
 * @param qmmmReferenceModels A vector that stores the QM/MM reference models.
 * @param fullStructure The molecular structure of the full system.
 * @param bondOrders Bond orders of the full system.
 * @param settings The settings.
 * @param log The logger.
 */
void generateQmRegionCandidates(std::vector<QmmmModel>& qmmmModelCandidates, std::vector<QmmmModel>& qmmmReferenceModels,
                                const Utils::AtomCollection& fullStructure, const Utils::BondOrderCollection& bondOrders,
                                const Utils::Settings& settings, Core::Log& log);

/**
 * @brief Given a cutting probability of 100 percent, one only gets a single QM region candidate
 *        for a given initial radius.
 * @param qmmmModelCandidates A vector of in which the QM/MM candidate models are stored.
 * @param generator An already prepared instance of a subsystem generator class.
 * @param fragmentAnalyzer An instance of the fragment analyzer that was passed to the generator.
 * @param initialRadius The initial radius for the the QM region generation.
 * @param centerAtoms A list of the center atom(s) around which the QM region(s) will be constructed.
 * @return The generated qmmm model.
 */
QmmmModel generateSingleQmRegionNonRandomly(SwooseUtilities::SubsystemGenerator& generator,
                                            SwooseUtilities::FragmentAnalyzer& fragmentAnalyzer, const double& initialRadius,
                                            std::vector<int> centerAtoms, const Utils::AtomCollection fullStructure,
                                            const std::vector<std::list<int>> listsOfNeighbors);

} // namespace QmRegionCandidateGenerator
} // namespace Qmmm
} // namespace Scine

#endif // SWOOSE_QMMM_QMREGIONCANDIDATEGENERATOR_H
