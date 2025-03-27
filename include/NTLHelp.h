#pragma once

#include <vector>
#include <optional>

#include "NTL/GF2E.h"
#include "NTL/vec_GF2.h"
#include "NTL/mat_GF2.h"

namespace Math {
    size_t GF2EToIndex(const NTL::GF2E& x);
    size_t GF2XToIndex(const NTL::GF2X& x);
    size_t GF2VecToIndex(const NTL::vec_GF2& vec);
    
    NTL::GF2E IndexToGF2E(size_t i, size_t deg);
    NTL::GF2X IndexToGF2X(size_t i, size_t deg);
    NTL::vec_GF2 IndexToGF2Vec(size_t i, size_t deg);

    std::optional<size_t> HighestOnePos(const NTL::vec_GF2& vec);

    NTL::mat_GF2 InvRectMatrix(const NTL::mat_GF2& mat);
    size_t RunDoubleGaussianElimination(NTL::mat_GF2& mat);

    NTL::mat_GF2 BuildArikanKernel(size_t numLayers);
    NTL::mat_GF2 TransformFreezingMatrix_Perm(
        const NTL::mat_GF2& freezingMatrix,
        const NTL::mat_GF2& transformMatrix
    );
}
