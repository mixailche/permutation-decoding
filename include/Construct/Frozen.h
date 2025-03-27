#pragma once

#include <vector>

namespace Construct {
    std::vector<bool> BuildFrozenSet(const std::vector<double>& errorProbs, size_t dimension);
    std::vector<bool> BuildFrozenSet_Bhattacharyya(double p, size_t nLayers, size_t dimension);
    std::vector<bool> BuildFrozenSet_GA(double snr, size_t nLayers, size_t dimension);
    std::vector<bool> BuildFrozenSet_BlockPermInvariant(
        size_t numLayer, size_t dimension, size_t blockSize,
        const std::vector<double>& errorProbs
    );
    std::vector<bool> BuildFrozenSet_AutomorphismInvariant(
        size_t numLayers, size_t dimension,
        const std::vector<size_t>& blockSizes,
        const std::vector<double>& errorProbs,
        size_t minWeight = 0
    );

    std::vector<double> BhattacharyyaParameters(double p, size_t nLayers);
    std::vector<double> DensityEvolutionGA(double snr, size_t nLayers);
} // namespace construct
