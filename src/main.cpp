#include <iostream>

#include "Construct/PolarSubcode.h"
#include "Construct/Perm.hpp"
#include "Construct/Frozen.h"
#include "Simulate.hpp"

#include "NTLHelp.h"
#include <unordered_set>

static void PrintErrorRate(size_t numSteps, size_t numErrors)
{
    if (numSteps % 100 == 0 && numSteps > 0) {
        std::cout << "steps=" << numSteps
            << ", errors=" << numErrors
            << ", FER=" << static_cast<double>(numErrors) / numSteps
            << "\n";
    }
}

int main(void)
{
    double snr;
    size_t length, dimension, numDFS_A, numDFS_B, maxPaths, permsCount;
    std::cin >> snr >> length >> dimension >> numDFS_A >> numDFS_B >> maxPaths >> permsCount;

    auto numLayers = Utils::IntLog2(length);

    constexpr size_t lastBlockSize = 3;
    std::vector<size_t> blockSizes;
    for (size_t i = 0; i < Utils::IntLog2(length) - lastBlockSize; i++) {
        blockSizes.push_back(1);
    }
    blockSizes.push_back(lastBlockSize);

    using namespace Construct;

    auto errorProbs = DensityEvolutionGA(2.0, numLayers);
    auto polar = BuildFrozenSet_AutomorphismInvariant(numLayers, dimension + numDFS_A, blockSizes, errorProbs);
    auto linearPerms = Utils::ListBlockDigitsPermutations(blockSizes);
    auto spec = BuildPermFriendlyRandomizedPolarSubcode(polar, errorProbs, linearPerms, numDFS_A, numDFS_B);
    auto perms = BuildBLTAPermSet(permsCount, numLayers, 1, 1, blockSizes);

    constexpr size_t numIterations = 10'000'000;
    constexpr size_t maxErrors = 100;

    Running::Simulate_PolarSubcodePermSCL(snr, numIterations, maxErrors, spec, perms, maxPaths, PrintErrorRate);

    return 0;
}
