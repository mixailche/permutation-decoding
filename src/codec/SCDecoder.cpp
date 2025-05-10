#include <algorithm>

#include "codec/SCDecoder.h"
#include "utils/Utils.hpp"

using codec::SCDecoder;

SCDecoder::SCDecoder(const PolarSpecification* spec)
    : mSpec(spec)
    , mNumLayers(utils::IntLog2(mSpec->Length))
    , mMetric(0)
    , mLLRs(mNumLayers + 1)
    , mSymbols(mNumLayers + 1)
    , mDynFrozenSymbols(mSpec->Length)
{
    for (size_t layer = 0; layer <= mNumLayers; layer++) {
        auto layerSize = GetLayerSize(layer);
        mLLRs[layer] = new double[layerSize];
        mSymbols[layer] = new bool[layerSize];
    }
}

SCDecoder::~SCDecoder()
{
    for (size_t layer = 0; layer <= mNumLayers; layer++) {
        delete[] mLLRs[layer];
        delete[] mSymbols[layer];
    }
}

double SCDecoder::Metric() const
{
    return mMetric;
}

size_t SCDecoder::GetLayerSize(size_t layer) const
{
    return 1ull << layer;
}

size_t SCDecoder::NumOperations() const
{
    return mNumOperations;
}

std::vector<bool> SCDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    mMetric = 0;
    mNumOperations = 0;
    std::fill(mDynFrozenSymbols.begin(), mDynFrozenSymbols.end(), 0);
    std::copy(inputLLRs.begin(), inputLLRs.end(), mLLRs[mNumLayers]);
    ProcessLayer(mNumLayers, 0);
    return std::vector(mSymbols[mNumLayers], mSymbols[mNumLayers] + mSpec->Length);
}

void SCDecoder::ProcessLayer(size_t layer, size_t phase) const
{
    if (layer == 0) {
        return ProcessPhase(phase);
    }

    auto half = GetLayerSize(layer - 1);

    // Negative branch
    for (size_t i = 0; i < half; i++) {
        auto a = mLLRs[layer][i];
        auto b = mLLRs[layer][i + half];
        mLLRs[layer - 1][i] = utils::NegativeMerge(a, b);
        mNumOperations++;
    }
    ProcessLayer(layer - 1, 2 * phase);

    // Positive branch
    for (size_t i = 0; i < half; i++) {
        auto a = mLLRs[layer][i];
        auto b = mLLRs[layer][i + half];
        auto u = mSymbols[layer - 1][i];
        mSymbols[layer][i] = u;
        mLLRs[layer - 1][i] = utils::PositiveMerge(a, b, u);
        mNumOperations++;
    }
    ProcessLayer(layer - 1, 2 * phase + 1);

    for (size_t i = 0; i < half; i++) {
        auto u = mSymbols[layer - 1][i];
        mSymbols[layer][i] ^= u;
        mSymbols[layer][i + half] = u;
    }
}

void SCDecoder::ProcessPhase(size_t phase) const
{
    auto y = mLLRs[0][0];

    if (mSpec->StaticFrozen[phase]) {
        AssignSymbol(phase, 0);
        mMetric += std::min(0.0, y);
        mNumOperations++;
    }
    else if (mSpec->Dynamic.Frozen[phase]) {
        auto symbol = mDynFrozenSymbols[phase];
        AssignSymbol(phase, symbol);
        if (symbol != y < 0) {
            mMetric -= std::abs(y);
            mNumOperations++;
        }
    }
    else {
        AssignSymbol(phase, y < 0);
    }
}

void SCDecoder::AssignSymbol(size_t phase, bool value) const
{
    mSymbols[0][0] = value;
    for (size_t dyn : mSpec->Dynamic.ForwardEquations[phase]) {
        mDynFrozenSymbols[dyn] = mDynFrozenSymbols[dyn] ^ value;
    }
}
