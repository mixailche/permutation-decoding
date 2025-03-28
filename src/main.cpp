#include <iostream>

#include "Construct/PolarSubcode.h"
#include "Construct/Perm.hpp"
#include "Construct/Frozen.h"
#include "Simulate.hpp"

#include "NTLHelp.h"

static void PrintErrorRate(size_t numSteps, size_t numErrors)
{
    if (numSteps % 1000 == 0 && numSteps > 0) {
        std::cout << "steps=" << numSteps
            << ", errors=" << numErrors
            << ", FER=" << static_cast<double>(numErrors) / numSteps;
    }
}

int main(void)
{
    // Dummy example
    double snr;
    size_t length, dimension, numDFS_A, numDFS_B, l, h, maxPaths, permsCount;
    std::cin >> snr >> length >> dimension >> numDFS_A >> numDFS_B >> l >> h >> maxPaths >> permsCount;
    
    auto numLayers = Utils::IntLog2(length);

    std::vector<size_t> blockSizesL = { l };
    for (size_t i = 0; i < Utils::IntLog2(numLayers) - l; i++) {
        blockSizesL.push_back(1);
    }

    std::vector<size_t> blockSizesH = { h };
    for (size_t i = 0; i < Utils::IntLog2(length) - h; i++) {
        blockSizesH.push_back(1);
    }

    using namespace Construct;

    auto errorProbs = DensityEvolutionGA(2.0, numLayers);
    auto polar = BuildFrozenSet_AutomorphismInvariant(numLayers, dimension, blockSizesH, errorProbs);
    auto linearPerms = Utils::ListBlockDigitsPermutations(blockSizesL);
    auto spec = BuildPermFriendlyRandomizedPolarSubcode(polar, errorProbs, linearPerms, numDFS_A, numDFS_B, l);
    auto perms = BuildJointPermSet(permsCount, numLayers, l, h);

    Running::Simulate_PolarSubcodePermSCL(snr, 10'000'000, 100, spec, perms, maxPaths, PrintErrorRate);

    return 0;
}
