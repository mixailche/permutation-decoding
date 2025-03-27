#include <algorithm>
#include "PolarCode/SCDecoder.h"
#include "Utils.hpp"

using Codec::PolarCode::SCDecoder;

SCDecoder::SCDecoder(const PolarCodeSpecification* spec)
    : mSpec(spec)
    , mNLayers(Utils::IntLog2(mSpec->Length))
    , mMetric(0)
    , mLLRs(mNLayers + 1)
    , mSymbols(mNLayers + 1)
{
    for (size_t layer = 0; layer <= mNLayers; layer++) {
        auto layerSize = GetLayerSize(layer);
        mLLRs[layer] = new double[layerSize];
        mSymbols[layer] = new bool[layerSize];
    }
}

SCDecoder::~SCDecoder()
{
    for (size_t layer = 0; layer <= mNLayers; layer++) {
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

Codec::InfVector SCDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    mMetric = 0;
    std::copy(inputLLRs.begin(), inputLLRs.end(), mLLRs[mNLayers]);
    ProcessLayer(mNLayers, 0);
    return std::vector(mSymbols[mNLayers], mSymbols[mNLayers] + mSpec->Length);
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
        mLLRs[layer - 1][i] = Utils::NegativeMerge(a, b);
    }
    ProcessLayer(layer - 1, 2 * phase);

    // Positive branch
    for (size_t i = 0; i < half; i++) {
        auto a = mLLRs[layer][i];
        auto b = mLLRs[layer][i + half];
        auto u = mSymbols[layer - 1][i];
        mSymbols[layer][i] = u;
        mLLRs[layer - 1][i] = Utils::PositiveMerge(a, b, u);
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

    if (mSpec->Frozen[phase]) {
        mSymbols[0][0] = 0;
        mMetric += std::min(0.0, y);
    }
    else {
        mSymbols[0][0] = y < 0;
    }
}
