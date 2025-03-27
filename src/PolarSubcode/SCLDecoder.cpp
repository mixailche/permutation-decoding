#include <numeric>
#include <algorithm>

#include "PolarSubcode/SCLDecoder.h"
#include "Utils.hpp"

using Codec::PolarSubcode::SCLDecoder;

SCLDecoder::SCLDecoder(const PolarSubcodeSpecification* spec, size_t maxPaths)
    : mSpec(spec)
    , mNLayers(Utils::IntLog2(spec->Length))
    , mMaxPaths(maxPaths)
    , mPathMetrics(maxPaths)
    , mPathLLRs(mNLayers + 1, std::vector<double*>(maxPaths))
    , mPathSymbols(mNLayers + 1, std::vector<bool*>(maxPaths))
    , mInactiveArrayIndices(mNLayers + 1)
    , mInactivePathIndices()
    , mIsActivePath(maxPaths, false)
    , mArrayReferenceCount(mNLayers + 1, std::vector<size_t>(maxPaths, 0))
    , mPathToArray(mNLayers + 1, std::vector<size_t>(maxPaths))
    , mPathDynFrozenSymbols(maxPaths)
{
    for (size_t layer = 0; layer <= mNLayers; layer++) {
        auto layerSize = GetLayerSize(layer);

        for (size_t pathIdx = 0; pathIdx < mMaxPaths; pathIdx++) {
            mPathLLRs[layer][pathIdx] = new double[layerSize];
            mPathSymbols[layer][pathIdx] = new bool[layerSize];
        }
    }

    for (size_t pathIdx = 0; pathIdx < maxPaths; pathIdx++) {
        mPathDynFrozenSymbols[pathIdx] = new bool[mSpec->Length];
    }
}

double SCLDecoder::Metric() const
{
    return -mBestMetric;
}

template <typename Container>
static inline void DeleteArrays(Container& c)
{
    std::for_each(c.begin(), c.end(), [](auto* arr) { delete[] arr; });
}

SCLDecoder::~SCLDecoder()
{
    DeleteArrays(mPathDynFrozenSymbols);

    for (size_t layer = 0; layer <= mNLayers; layer++) {
        DeleteArrays(mPathLLRs[layer]);
        DeleteArrays(mPathSymbols[layer]);
    }
}

inline size_t SCLDecoder::GetLayerSize(size_t layer) const
{
    return 1ull << layer;
}

void SCLDecoder::Initialize() const
{
    std::fill(mPathMetrics.begin(), mPathMetrics.end(), 0);
    std::fill(mIsActivePath.begin(), mIsActivePath.end(), false);

    for (auto& layerReferenceCount : mArrayReferenceCount) {
        std::fill(layerReferenceCount.begin(), layerReferenceCount.end(), 0);
    }

    for (auto dynFrozenSymbols : mPathDynFrozenSymbols) {
        std::fill_n(dynFrozenSymbols, mSpec->Length, 0);
    }

    for (size_t layer = 0; layer <= mNLayers; layer++) {
        mInactiveArrayIndices[layer].clear();
        for (size_t pathIdx = 0; pathIdx < mMaxPaths; pathIdx++) {
            mInactiveArrayIndices[layer].push_back(pathIdx);
        }
    }

    mInactivePathIndices.clear();
    for (size_t pathIdx = 0; pathIdx < mMaxPaths - 1; pathIdx++) {
        mInactivePathIndices.push_back(pathIdx);
    }

    auto activePathIdx = mMaxPaths - 1;
    mIsActivePath[activePathIdx] = true;

    for (size_t layer = 0; layer <= mNLayers; mInactiveArrayIndices[layer++].pop_back()) {
        auto arrayIdx = mInactiveArrayIndices[layer].back();
        mPathToArray[layer][activePathIdx] = arrayIdx;
        mArrayReferenceCount[layer][arrayIdx] = 1;
    }
}

size_t SCLDecoder::ClonePath(size_t pathIdx) const
{
    auto newPathIdx = mInactivePathIndices.back(); mInactivePathIndices.pop_back();
    mPathMetrics[newPathIdx] = mPathMetrics[pathIdx];
    mIsActivePath[newPathIdx] = true;

    for (size_t layer = 0; layer <= mNLayers; layer++) {
        auto arrayIdx = mPathToArray[layer][pathIdx];
        mPathToArray[layer][newPathIdx] = arrayIdx;
        mArrayReferenceCount[layer][arrayIdx]++;
    }

    std::copy_n(mPathDynFrozenSymbols[pathIdx], mSpec->Length, mPathDynFrozenSymbols[newPathIdx]);

    return newPathIdx;
}

void SCLDecoder::KillPath(size_t pathIdx) const
{
    mIsActivePath[pathIdx] = false;
    mInactivePathIndices.push_back(pathIdx);

    for (size_t layer = 0; layer <= mNLayers; layer++) {
        auto array_idx = mPathToArray[layer][pathIdx];
        if (--mArrayReferenceCount[layer][array_idx] == 0) {
            mInactiveArrayIndices[layer].push_back(array_idx);
        }
    }
}

size_t SCLDecoder::GetArrayIndex(size_t layer, size_t pathIdx) const
{
    auto srcIdx = mPathToArray[layer][pathIdx];

    if (mArrayReferenceCount[layer][srcIdx] == 1) {
        return srcIdx;
    }

    auto layerSize = GetLayerSize(layer);
    auto dstIdx = mInactiveArrayIndices[layer].back();
    mInactiveArrayIndices[layer].pop_back();

    std::copy_n(mPathLLRs[layer][srcIdx], layerSize, mPathLLRs[layer][dstIdx]);
    std::copy_n(mPathSymbols[layer][srcIdx], layerSize, mPathSymbols[layer][dstIdx]);

    mArrayReferenceCount[layer][srcIdx]--;
    mArrayReferenceCount[layer][dstIdx] = 1;

    mPathToArray[layer][pathIdx] = dstIdx;

    return dstIdx;
}

double* SCLDecoder::GetLLRsArray(size_t layer, size_t pathIdx) const
{
    auto arrayIdx = GetArrayIndex(layer, pathIdx);
    return mPathLLRs[layer][arrayIdx];
}

bool* SCLDecoder::GetSymbolsArray(size_t layer, size_t pathIdx) const
{
    auto arrayIdx = GetArrayIndex(layer, pathIdx);
    return mPathSymbols[layer][arrayIdx];
}

void SCLDecoder::ProcessLayer(size_t layer, size_t phase) const
{
    if (layer == 0) {
        return ProcessPhase(phase);
    }

    ProcessNegativeBranch(layer, phase);
    ProcessPositiveBranch(layer, phase);

    for (size_t pathIdx = 0; pathIdx < mMaxPaths; pathIdx++) {
        if (!mIsActivePath[pathIdx]) {
            continue;
        }

        auto currentSymbols = GetSymbolsArray(layer, pathIdx);
        auto prevSymbols = GetSymbolsArray(layer - 1, pathIdx);

        for (size_t i = 0, half = GetLayerSize(layer - 1); i < half; i++) {
            auto u = prevSymbols[i];
            currentSymbols[i] ^= u;
            currentSymbols[i + half] = u;
        }
    }
}

void SCLDecoder::ProcessNegativeBranch(size_t layer, size_t phase) const
{
    for (size_t pathIdx = 0; pathIdx < mMaxPaths; pathIdx++) {
        if (!mIsActivePath[pathIdx]) {
            continue;
        }

        auto currentLLRs = GetLLRsArray(layer, pathIdx);
        auto prevLLRs = GetLLRsArray(layer - 1, pathIdx);

        auto currentSymbols = GetSymbolsArray(layer, pathIdx);
        auto prevSymbols = GetSymbolsArray(layer - 1, pathIdx);

        for (size_t i = 0, half = GetLayerSize(layer - 1); i < half; i++) {
            auto a = currentLLRs[i];
            auto b = currentLLRs[i + half];
            prevLLRs[i] = Utils::NegativeMerge(a, b);
        }
    }

    ProcessLayer(layer - 1, 2 * phase);
}

void SCLDecoder::ProcessPositiveBranch(size_t layer, size_t phase) const
{
    for (size_t pathIdx = 0; pathIdx < mMaxPaths; pathIdx++) {
        if (!mIsActivePath[pathIdx]) {
            continue;
        }

        auto currentLLRs = GetLLRsArray(layer, pathIdx);
        auto prevLLRs = GetLLRsArray(layer - 1, pathIdx);

        auto currentSymbols = GetSymbolsArray(layer, pathIdx);
        auto prevSymbols = GetSymbolsArray(layer - 1, pathIdx);

        for (size_t i = 0, half = GetLayerSize(layer - 1); i < half; i++) {
            auto a = currentLLRs[i];
            auto b = currentLLRs[i + half];
            auto u = prevSymbols[i];
            currentSymbols[i] = u;
            prevLLRs[i] = Utils::PositiveMerge(a, b, u);
        }
    }

    ProcessLayer(layer - 1, 2 * phase + 1);
}

void SCLDecoder::ProcessPhase(size_t phase) const
{
    if (mSpec->StaticFrozen[phase]) {
        ContinuePaths_StaticFrozen(phase);
    }
    else if (mSpec->Dynamic.Frozen[phase]) {
        ContinuePaths_DynamicFrozen(phase);
    }
    else {
        ContinuePaths_Unfrozen(phase);
    }
}

static double GetContinuedPathMetric(double currentMetric, double y, bool u)
{
    return u == y < 0 ? currentMetric : currentMetric + std::abs(y);
}

void SCLDecoder::ContinuePaths_StaticFrozen(size_t phase) const
{
    for (size_t pathIdx = 0; pathIdx < mMaxPaths; pathIdx++) {
        if (mIsActivePath[pathIdx]) {
            AssignSymbol(pathIdx, phase, 0);
        }
    }
}

void SCLDecoder::ContinuePaths_DynamicFrozen(size_t phase) const
{
    for (size_t pathIdx = 0; pathIdx < mMaxPaths; pathIdx++) {
        if (mIsActivePath[pathIdx]) {
            AssignSymbol(pathIdx, phase, mPathDynFrozenSymbols[pathIdx][phase]);
        }
    }
}

void SCLDecoder::ContinuePaths_Unfrozen(size_t phase) const
{
    std::vector<double> forks(2 * mMaxPaths);

    for (size_t i = 0; i < mMaxPaths; i++) {
        if (!mIsActivePath[i]) {
            forks[2 * i] = forks[2 * i + 1] = std::numeric_limits<double>::max();
        }
        else {
            auto pathLLRs = GetLLRsArray(0, i);
            forks[2 * i] = GetContinuedPathMetric(mPathMetrics[i], pathLLRs[0], 0);
            forks[2 * i + 1] = GetContinuedPathMetric(mPathMetrics[i], pathLLRs[0], 1);
        }
    }

    size_t nActivePaths = std::count(mIsActivePath.begin(), mIsActivePath.end(), true);
    size_t nContinuedPaths = std::min(mMaxPaths, 2 * nActivePaths);

    auto indices = Utils::SortingPerm(forks, std::less<double>());
    std::vector<Utils::pair<bool>> continuations(mMaxPaths, { false, false });

    for (size_t i = 0; i < nContinuedPaths; i++) {
        auto pathIdx = indices[i] / 2;
        auto u = indices[i] % 2;
        continuations[pathIdx][u] = true;
    }

    for (size_t pathIdx = 0; pathIdx < mMaxPaths; pathIdx++) {
        auto& cont = continuations[pathIdx];
        if (mIsActivePath[pathIdx] && !cont[0] && !cont[1]) {
            KillPath(pathIdx);
        }
    }

    for (size_t pathIdx = 0; pathIdx < mMaxPaths; pathIdx++) {
        auto& cont = continuations[pathIdx];

        if (!cont[0] && !cont[1]) {
            continue;
        }

        if (cont[0] && cont[1]) {
            auto newPathIdx = ClonePath(pathIdx);
            AssignSymbol(newPathIdx, phase, 0);
            AssignSymbol(pathIdx, phase, 1);
        }
        else if (cont[0]) {
            AssignSymbol(pathIdx, phase, 0);
        }
        else {
            AssignSymbol(pathIdx, phase, 1);
        }
    }
}

void SCLDecoder::AssignSymbol(size_t pathIdx, size_t phase, bool symbol) const
{
    auto pathLLRs = GetLLRsArray(0, pathIdx);
    auto pathSymbols = GetSymbolsArray(0, pathIdx);
    auto pathMetric = mPathMetrics[pathIdx];

    pathSymbols[0] = symbol;
    mPathMetrics[pathIdx] = GetContinuedPathMetric(pathMetric, pathLLRs[0], symbol);

    for (auto dyn : mSpec->Dynamic.ForwardEquations[phase]) {
        mPathDynFrozenSymbols[pathIdx][dyn] ^= symbol;
    }
}

Codec::InfVector SCLDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    Initialize();

    auto activePathLLRs = GetLLRsArray(mNLayers, mMaxPaths - 1);
    std::copy(inputLLRs.begin(), inputLLRs.end(), activePathLLRs);

    ProcessLayer(mNLayers, 0);

    size_t bestPathIdx = -1;
    mBestMetric = std::numeric_limits<double>::max();

    for (size_t pathIdx = 0; pathIdx < mMaxPaths; pathIdx++) {
        auto pathMetric = mPathMetrics[pathIdx];
        if (mIsActivePath[pathIdx] && pathMetric < mBestMetric) {
            bestPathIdx = pathIdx;
            mBestMetric = pathMetric;
        }
    }

    auto bestPathSymbols = GetSymbolsArray(mNLayers, bestPathIdx);
    return std::vector(bestPathSymbols, bestPathSymbols + mSpec->Length);
}
