#pragma once

#include <ostream>
#include <vector>
#include <optional>

#include "GField.h"
#include "VecGF2.h"
#include "MatGF2.h"

namespace math {
    size_t VecGF2ToIndex(const VecGF2& vec);
    VecGF2 IndexToVecGF2(size_t index, size_t dimension);
    
    size_t HammingDistance(const VecGF2& lhs, const VecGF2& rhs);

    MatGF2 RandomMatGF2(size_t numRows, size_t numCols);
    VecGF2 RandomVecGF2(size_t dimension);

    std::ostream& operator<<(std::ostream& ostr, const MatGF2& mat);
    std::ostream& operator<<(std::ostream& ostr, const VecGF2& vec);

    MatGF2 BuildArikanKernel(size_t numLayers);
    
    MatGF2 TransformFreezingMatrix(const MatGF2& freezingMatrix, const MatGF2& permMatrix);
}
