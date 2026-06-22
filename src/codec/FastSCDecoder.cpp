#include <algorithm>

#include "codec/FastSCDecoder.h"
#include "utils/Utils.hpp"

using codec::FastSCDecoder;

FastSCDecoder::FastSCDecoder(const PolarSpecification* spec)
    : mSpec(spec)
    , mNumLayers(utils::IntLog2(mSpec->Length))
    , mMetric(0)
    , mLLRs(mNumLayers + 1)
    , mSymbols(mNumLayers + 1)
    , mDynFrozenSymbols(mSpec->Length)
    , mCodeword(mSpec->Length)
{
    for (size_t layer = 0; layer <= mNumLayers; layer++) {
        auto layerSize = GetLayerSize(layer);
        mLLRs[layer] = new double[layerSize];
        mSymbols[layer] = new bool[layerSize];
    }
}

FastSCDecoder::~FastSCDecoder()
{
    for (size_t layer = 0; layer <= mNumLayers; layer++) {
        delete[] mLLRs[layer];
        delete[] mSymbols[layer];
    }
}

double FastSCDecoder::Metric() const
{
    return mMetric;
}

size_t FastSCDecoder::GetLayerSize(size_t layer) const
{
    return 1ull << layer;
}

size_t FastSCDecoder::NumOperations() const
{
    return mNumOperations;
}

std::vector<bool> FastSCDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    mMetric = 0;
    mNumOperations = 0;
    std::fill(mDynFrozenSymbols.begin(), mDynFrozenSymbols.end(), 0);
    std::copy(inputLLRs.begin(), inputLLRs.end(), mLLRs[mNumLayers]);
    ProcessLayer(mNumLayers, 0);
    return std::vector(mSymbols[mNumLayers], mSymbols[mNumLayers] + mSpec->Length);
}

static void Polarize(bool* begin, bool* end)
{
    auto size = std::distance(begin, end);
    if (size == 1) {
        return;
    }

    auto half = size / 2;
    auto middle = begin + half;

    Polarize(begin, middle);
    Polarize(middle, end);

    std::transform(begin, middle, middle, begin, [](bool lhs, bool rhs) {
        return lhs ^ rhs;
    });
}

void FastSCDecoder::ProcessLayer(size_t layer, size_t phase) const
{
    auto length = GetLayerSize(layer);
    auto half = GetLayerSize(layer - 1);

    for (size_t i = 0; i < length; i++) {
        mSymbols[layer][i] = mLLRs[layer][i] < 0;
    }
    Polarize(mSymbols[layer], mSymbols[layer] + length);
    
    bool parityCheck = true;
    for (size_t i = 0; i < length; i++) {
        if (mSpec->StaticFrozen[phase * (1ull << layer) + i] && mSymbols[layer][i] == 1) {
            parityCheck = false;
            break;
        }
    }
    if (parityCheck) {
        Polarize(mSymbols[layer], mSymbols[layer] + length);
        return;
    }

    if (layer == 0) {
        return ProcessPhase(phase);
    }

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

void FastSCDecoder::ProcessPhase(size_t phase) const
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

void FastSCDecoder::AssignSymbol(size_t phase, bool value) const
{
    mSymbols[0][0] = value;
    for (size_t dyn : mSpec->Dynamic.ForwardEquations[phase]) {
        mDynFrozenSymbols[dyn] = mDynFrozenSymbols[dyn] ^ value;
    }
}

//bool FastSCDecoder::CheckStopCondition(size_t layer, size_t phase) const
//{
//    if (layer == mNumLayers) {
//        // ...
//        return;
//    }
//
//    auto length = 1ull << layer;
//
//    if (phase % 2 == 0) {
//        for (size_t i = 0; i < length; i++) {
//            mNumOperations++;
//            auto rightLLR = utils::PositiveMerge(mLLRs[layer - 1][i], mLLRs[layer - 1][i + length], mCodeword[i]);
//            bool rightSymbol = rightLLR < 0 ? 1 : 0;
//            mCodeword[i] = mCodeword[i] ^ rightSymbol;
//            mCodeword[i + length] = rightSymbol;
//        }
//    }
//    else {
//        for (size_t i = 0; i < length; i++) {
//            mCodeword[i + length] = mCodeword[i];
//            mCodeword[i] = mSymbols[layer][i] ^ mCodeword[i];
//        }
//    }
//}
