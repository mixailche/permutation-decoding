#pragma once

#include "math/MatGF2.h"
#include "math/VecGF2.h"
#include "math/GField.h"
#include "math/Perm.hpp"
#include "utils/Utils.hpp"
#include "codec/PolarSpecification.h"

namespace construct {
    codec::PolarSpecification BuildEBCHSubcode(
        size_t numLayers, size_t dimension, size_t minDist,
        const math::GField& field,
        const std::vector<double>& errorProbs
    );

    codec::PolarSpecification BuildPermFriendlyEBCHSubcode(
        size_t numLayers, size_t dimension, size_t minDist,
        const math::GField& field,
        const std::vector<double>& errorProbs,
        const std::vector<size_t>& blockSizes
    );

    codec::PolarSpecification BuildRandomizedPolarSubcode(
        const std::vector<bool>& baseCodeFrozenSymbols,
        const std::vector<double>& errorProbs, size_t numDFS_A, size_t numDFS_B
    );

    codec::PolarSpecification BuildPermFriendlyRandomizedPolarSubcode(
        const std::vector<bool>& baseCodeFrozenSymbols,
        const std::vector<double>& errorProbs,
        const std::vector<utils::RawPerm>& linearPerms,
        size_t numDFS_A, size_t numDFS_B
    );

    codec::PolarSpecification BuildPlotkinPolarSubcode(
        const codec::PolarSpecification& code1,
        const codec::PolarSpecification& code2
    );
}
