#pragma once

#include <cstddef>
#include <vector>

#include "NTL/GF2X.h"
#include "NTL/mat_GF2.h"
#include "Specification.h"
#include "Perm.hpp"
#include "Utils.hpp"

namespace Construct {
    NTL::mat_GF2 BuildFreezingMatrix(const Codec::PolarSubcodeSpecification& spec);

    std::vector<std::vector<size_t>> BuildDynamicFreezingEquations(
        const NTL::mat_GF2& freezingMatrix,
        std::vector<bool>& staticFrozen
    );

    Codec::PolarSubcodeSpecification BuildEBCHSubcode(
        size_t numLayers, size_t dimension, size_t minDist,
        const std::vector<double>& errorProbs
    );

    Codec::PolarSubcodeSpecification BuildAutFriendlyEBCHSubcode(
        size_t numLayers, size_t dimension, size_t minDist,
        const std::vector<double>& errorProbs,
        const std::vector<size_t>& blockSizes
    );

    Codec::PolarSubcodeSpecification BuildCyclicPolarSubcode(
        size_t numLayers, size_t dimension,
        const std::vector<double>& errorProbs
    );

    Codec::PolarSubcodeSpecification BuildRandomizedPolarSubcode(
        const Codec::PolarCodeSpecification& polarSpec,
        const std::vector<double>& errorProbs, size_t numDFS_A, size_t numDFS_B
    );

    Codec::PolarSubcodeSpecification BuildPermFriendlyRandomizedPolarSubcode(
        const Codec::PolarCodeSpecification& polarSpec,
        const std::vector<double>& errorProbs,
        const std::vector<Utils::RawPerm>& linearPerms,
        size_t numDFS_A, size_t numDFS_B, size_t l
    );

    Codec::PolarSubcodeSpecification BuildPlotkinPolarSubcode(
        const Codec::PolarSubcodeSpecification& code1,
        const Codec::PolarSubcodeSpecification& code2
    );
}
