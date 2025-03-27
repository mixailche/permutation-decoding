#include <vector>
#include <numeric>
#include <queue>
#include <unordered_map>

#include "Utils.hpp"
#include "NTLHelp.h"

size_t Math::GF2EToIndex(const NTL::GF2E& x)
{
    return GF2XToIndex(NTL::conv<NTL::GF2X>(x));
}

size_t Math::GF2XToIndex(const NTL::GF2X& x)
{
    return GF2VecToIndex(NTL::to_vec_GF2(x));
}

size_t Math::GF2VecToIndex(const NTL::vec_GF2& vec)
{
    size_t i = 0;
    for (size_t j = 0; j < vec.length(); j++) {
        if (NTL::IsOne(vec[j])) {
            i += 1ull << j;
        }
    }
    return i;
}

NTL::GF2E Math::IndexToGF2E(size_t i, size_t deg)
{
    return NTL::to_GF2E(IndexToGF2X(i, deg));
}

NTL::GF2X Math::IndexToGF2X(size_t i, size_t deg)
{
    return NTL::to_GF2X(IndexToGF2Vec(i, deg));
}

NTL::vec_GF2 Math::IndexToGF2Vec(size_t i, size_t deg)
{
    NTL::vec_GF2 vec(NTL::INIT_SIZE, deg);
    for (size_t j = 0; j < deg; j++) {
        if (i & (1ull << j)) {
            vec[j] = 1;
        }
    }
    return vec;
}

static bool IsIndependentSet(
    const NTL::vec_vec_GF2& vectors, size_t numVectors, size_t length)
{
    auto mat = NTL::ident_mat_GF2(length);
    for (size_t i = 0; i < numVectors; i++) {
        mat[i] = vectors[i];
    }
    return NTL::IsOne(NTL::determinant(mat));
}

static inline size_t LowestOnePos(const NTL::vec_GF2& vec)
{
   for (size_t i = 0; i < vec.length(); i++) {
        if (NTL::IsOne(vec[i])) {
            return i;
        }
    }
    return vec.length();
}

NTL::mat_GF2 Math::InvRectMatrix(const NTL::mat_GF2& mat)
{
    // numRows < numCols
    size_t numRows = mat.NumRows();
    size_t numCols = mat.NumCols();

    NTL::mat_GF2 echelonForm = mat;
    NTL::gauss(echelonForm);

    std::vector<size_t> colsPerm(numCols);
    std::vector usedCols(numCols, false);

    for (size_t i = 0; i < numRows; i++) {
        colsPerm[i] = LowestOnePos(echelonForm[i]);
        usedCols[colsPerm[i]] = true;
    }

    for (size_t i = numRows, j = 0; i < numCols; i++) {
        while (usedCols[j]) {
            j++;
        }
        colsPerm[i] = j;
        usedCols[j] = true;
    }

    auto extendedMat = NTL::ident_mat_GF2(numCols);
    for (size_t i = 0; i < numRows; i++) {
        for (size_t j = 0; j < numCols; j++) {
            extendedMat[i][j] = mat[i][colsPerm[j]];
        }
    }

    NTL::mat_GF2 inv = NTL::inv(extendedMat);
    NTL::mat_GF2 result(NTL::INIT_SIZE, numCols, numRows);
    for (size_t i = 0; i < result.NumRows(); i++) {
        for (size_t j = 0; j < result.NumCols(); j++) {
            result[colsPerm[i]][j] = inv[i][j];
        }
    }

    return result;
}

std::optional<size_t> Math::HighestOnePos(const NTL::vec_GF2& vec)
{
    std::optional<size_t> result = std::nullopt;
    for (size_t i = 0; i < vec.length(); i++) {
        if (NTL::IsOne(vec[i])) {
            result = i;
        }
    }
    return result;
}

size_t Math::RunDoubleGaussianElimination(NTL::mat_GF2& mat)
{
    std::vector<std::queue<size_t>> buckets(mat.NumCols());

    for (size_t r = 0; r < mat.NumRows(); r++) {
        if (auto c = HighestOnePos(mat[r]); c.has_value()) {
            buckets[*c].push(r);
        }
    }

    // { r, HighestOnePos(mat[r]) }, sorted by second descending
    std::vector<std::pair<size_t, size_t>> rowEnds;

    for (size_t c = mat.NumCols(); c-- > 0;) {
        auto& bucket = buckets[c];
        if (bucket.empty()) {
            continue;
        }

        auto r = bucket.front(); bucket.pop();
        auto& base = mat[r];
        rowEnds.emplace_back(r, c);

        while (!bucket.empty()) {
            auto r1 = bucket.front(); bucket.pop();
            auto& row = mat[r1];
            row += base;
            if (auto c1 = HighestOnePos(row); c1.has_value()) {
                buckets[*c1].push(r1);
            }
        }
    }

    for (size_t i = rowEnds.size(); i-- > 0;) {
        auto& [r, c] = rowEnds[i];
        auto& row = mat[r];

        for (size_t j = 0; j < i; j++) {
            auto r1 = rowEnds[j].first;
            if (NTL::IsOne(mat[r1][c])) {
                mat[r1] += row;
            }
        }
    }

    return rowEnds.size();
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

NTL::mat_GF2 Math::BuildArikanKernel(size_t numLayers)
{
    static std::unordered_map<size_t, NTL::mat_GF2> cache;
    if (auto it = cache.find(numLayers); it != cache.end()) {
        return it->second;
    }

    auto length = 1ull << numLayers;
    auto elems = BuildArikanKernelElems(numLayers);
    auto result = NTL::ident_mat_GF2(length);

    for (size_t i = 0; i < length; i++) {
        for (size_t j = 0; j < i; j++) {
            result.put(i, j, elems[i][j]);
        }
    }

    cache.emplace(numLayers, result).second;
    return result;
}

NTL::mat_GF2 Math::TransformFreezingMatrix_Perm(
    const NTL::mat_GF2& freezingMatrix,
    const NTL::mat_GF2& permMatrix)
{
    auto numLayers = Utils::IntLog2(permMatrix.NumCols());

    //std::cout << transformMatrix.NumCols() << std::endl;
    auto kernel = NTL::transpose(BuildArikanKernel(numLayers));

    auto transformed = freezingMatrix * kernel * permMatrix * kernel;
    RunDoubleGaussianElimination(transformed);

    return transformed;
}
