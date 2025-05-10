#include <stdexcept>
#include <vector>
#include <random>
#include <unordered_map>

#include "utils/Utils.hpp"
#include "math/MathUtils.h"

using math::MatGF2;
using math::VecGF2;

size_t math::VecGF2ToIndex(const VecGF2& vec)
{
    //if (vec.Dimension() > 8 /* bits per byte */ * sizeof(size_t)) {
    //    throw std::invalid_argument("Vector is too long");
    //}

    size_t index = 0;
    for (size_t i = 0; i < vec.Dimension(); i++) {
        index |= vec[i] * (1ull << i);
    }
    return index;
}

VecGF2 math::IndexToVecGF2(size_t index, size_t dimension)
{
    std::vector<bool> vec(dimension);
    for (size_t i = 0; i < dimension; i++) {
        vec[i] = index & (1ull << i);
    }
    return vec;
}

size_t math::HammingDistance(const VecGF2& lhs, const VecGF2& rhs)
{
    return (lhs + rhs).HammingWeight();
}

MatGF2 math::RandomMatGF2(size_t numRows, size_t numCols)
{
    auto mat = MatGF2::Zeros(numRows, numCols);
    for (size_t i = 0; i < numRows; i++) {
        for (size_t j = 0; j < numCols; j++) {
            mat.Set(i, j, utils::RandomGF2());
        }
    }
    return mat;
}

VecGF2 math::RandomVecGF2(size_t dimension)
{
    VecGF2 vec(dimension, 0);
    for (size_t i = 0; i < dimension; i++) {
        vec.Set(i, utils::RandomGF2());
    }
    return vec;
}

std::ostream& math::operator<<(std::ostream& ostr, const math::MatGF2& mat)
{
    for (size_t i = 0; i < mat.NumRows(); i++) {
        ostr << mat[i];
        if (i != mat.NumRows()) {
            ostr << "\n";
        }
    }
    return ostr;
}

std::ostream& math::operator<<(std::ostream& ostr, const VecGF2& vec)
{
    ostr << "[ ";
    for (size_t i = 0; i < vec.Dimension(); i++) {
        ostr << vec[i] << " ";
    }
    return ostr << "]";
}

static std::vector<std::vector<bool>> BuildArikanKernelElems(size_t nLayers)
{
    if (nLayers == 0) {
        return { { 1 } };
    }

    auto size = 1ull << nLayers;
    auto half = size / 2;

    auto kernel = BuildArikanKernelElems(nLayers - 1);
    for (size_t i = 0; i < half; i++) {
        kernel.emplace_back(size);
        kernel[i].insert(kernel[i].end(), half, 0);

        for (size_t j = 0; j < half; j++) {
            auto symbol = kernel[i][j];
            kernel[half + i][j] = symbol;
            kernel[half + i][half + j] = symbol;
        }
    }

    return kernel;
}

MatGF2 math::BuildArikanKernel(size_t numLayers)
{
    static std::unordered_map<size_t, MatGF2> cache;
    if (auto it = cache.find(numLayers); it != cache.end()) {
        return it->second;
    }

    auto length = 1ull << numLayers;
    auto elems = BuildArikanKernelElems(numLayers);
    auto result = MatGF2::Identity(length);

    for (size_t i = 0; i < length; i++) {
        for (size_t j = 0; j < i; j++) {
            result.Set(i, j, elems[i][j]);
        }
    }

    cache.emplace(numLayers, result).second;
    return result;
}

MatGF2 math::TransformFreezingMatrix(const MatGF2& freezingMatrix, const MatGF2& permMatrix)
{
    auto numLayers = utils::IntLog2(permMatrix.NumCols());
    auto kernel = BuildArikanKernel(numLayers).Transpose();

    auto transformed = freezingMatrix * kernel * permMatrix * kernel;
    transformed.RunGaussianElimination();

    return transformed;
}
