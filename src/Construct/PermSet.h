#pragma once

#include <vector>
#include "math/Perm.hpp"
#include "codec/PolarSpecification.h"
#include "utils/Utils.hpp"

namespace construct {
    std::vector<size_t> GetStabBlocksStructure(const std::vector<bool>& frozen);

    std::vector<math::Perm> BuildDigitsPermSet_Random(
        size_t count, size_t numLayers, size_t numFrozenLayers, size_t minDist
    );

    std::vector<math::Perm> BuildDigitsPermSet_LeastErrorProb(
        size_t count, size_t numLayers, size_t minDist,
        const std::vector<double>& errorProbs,
        const std::vector<bool>& frozen
    );

    std::vector<math::Perm> BuildDigitsPermSet_LeastErrorProb(
        size_t count, size_t minDist,
        const std::vector<double>& errorProbs,
        const std::vector<size_t>& blockSizes,
        const codec::PolarSpecification& spec
    );

    std::vector<math::Perm> BuildJointPermSet(size_t count, size_t numLayers, size_t l, size_t h);

    std::vector<math::Perm> BuildBLTAPermSet(
        size_t count, size_t numLayers,
        size_t minDistUTElems, size_t minDistPerm,
        const std::vector<size_t>& blockSizes
    );
}
