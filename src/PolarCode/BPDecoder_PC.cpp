#include "Utils.hpp"
#include "PolarCode/BPDecoder.h"

using Codec::PolarCode::BPDecoder;

BPDecoder::BPDecoder(const PolarCodeSpecification* spec, size_t numIterMax)
    : BPDecoder(spec, numIterMax, Utils::Iota(Utils::IntLog2(spec->Length) + 1))
{}

BPDecoder::BPDecoder(const PolarCodeSpecification* spec, size_t numIterMax, std::vector<size_t> schedule)
    : mSpec(spec)
    , mNumLayers(Utils::IntLog2(mSpec->Length))
    , mNumIterMax(numIterMax)
    , mSchedule(std::move(schedule))
    , mEncoder(spec)
    , mLeftMessages(mNumLayers + 1)
    , mRightMessages(mNumLayers + 1)
{
    if (mNumIterMax == 0) {
        throw std::invalid_argument("Value of numIterMax must be > 0");
    }

    for (size_t layer = 0; layer <= mNumLayers; layer++) {
        mLeftMessages[layer] = new double[mSpec->Length];
        mRightMessages[layer] = new double[mSpec->Length];
    }
}

BPDecoder::~BPDecoder()
{
    for (size_t layer = 0; layer <= mNumLayers; layer++) {
        delete[] mLeftMessages[layer];
        delete[] mRightMessages[layer];
    }
}

inline size_t Codec::PolarCode::BPDecoder::GetLayerSize(size_t layer) const
{
    return 1ull << layer;
}

void BPDecoder::Initialize(const std::vector<double>& inputLLRs) const
{
    for (size_t layer = 0; layer <= mNumLayers; layer++) {
        std::fill_n(mLeftMessages[layer], mSpec->Length, 0);
        std::fill_n(mRightMessages[layer], mSpec->Length, 0);
    }

    for (size_t i = 0; i < mSpec->Length; i++) {
        mLeftMessages[mNumLayers][i] = inputLLRs[i];
        if (mSpec->Frozen[i]) {
            mRightMessages[mSchedule[0]][i] = std::numeric_limits<double>::max();
        }
    }
}

Codec::InfVector BPDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    Initialize(inputLLRs);
    Codec::InfVector codeword(inputLLRs.size());

    for (size_t i = 0; i < mNumIterMax; i++) {
        RunBPIteration();
        if (i % (mNumIterMax / 10) == 0 && CheckStopCondition(codeword)) {
            break;
        }
    }
    
    return GetCodeword();
}

std::pair<bool, std::vector<double>>
BPDecoder::SoftDecode(const std::vector<double>& inputLLRs) const
{
    Initialize(inputLLRs);
    Codec::InfVector codeword;

    for (size_t i = 0; i < mNumIterMax; i++) {
        RunBPIteration();
        if (i % (mNumIterMax / 10) == 0 && CheckStopCondition()) {
            return { true, GetLayerLLRs(mNumLayers) };
        }
    }

    return { false, GetLayerLLRs(mNumLayers) };
}

std::vector<std::pair<size_t, size_t>> BPDecoder::GetSymbolNodes(size_t layer) const
{
    auto layerSize = GetLayerSize(layer);
    std::vector used(mSpec->Length, false);
    std::vector<std::pair<size_t, size_t>> result;

    for (size_t i = 0; i < mSpec->Length; i++) {
        if (used[i]) {
            continue;
        }
        result.emplace_back(i, i + layerSize);
        used[i + layerSize] = true;
    }

    return result;
}

void BPDecoder::GetInputMessages(
    size_t leftLayer, size_t rightLayer, size_t i1, size_t i2,
    double& l1, double& l2, double& r1, double& r2
) const {
    l1 = mLeftMessages[rightLayer][i1];
    l2 = mLeftMessages[rightLayer][i2];
    r1 = mRightMessages[leftLayer][i1];
    r2 = mRightMessages[leftLayer][i2];
}

void BPDecoder::RunBPIteration() const
{
    double l1, l2, r1, r2;

    for (size_t i = mNumLayers; i-- > 0;) {
        auto currentLayer = mSchedule[i];
        auto prevLayer = (i == mNumLayers - 1) ? mNumLayers : mSchedule[i + 1];
        auto symbolNodes = GetSymbolNodes(currentLayer);
        for (auto& [i1, i2] : symbolNodes) {
            GetInputMessages(currentLayer, prevLayer, i1, i2, l1, l2, r1, r2);
            mLeftMessages[currentLayer][i1] = Utils::BoxPlus(l1, l2 + r2);
            mLeftMessages[currentLayer][i2] = Utils::BoxPlus(l1, r1) + l2;
        }
    }

    for (size_t i = 0; i < mNumLayers; i++) {
        auto currentLayer = (i == mNumLayers - 1) ? mNumLayers : mSchedule[i + 1];
        auto prevLayer = mSchedule[i];
        auto symbolNodes = GetSymbolNodes(prevLayer);
        for (auto& [i1, i2] : symbolNodes) {
            GetInputMessages(prevLayer, currentLayer, i1, i2, l1, l2, r1, r2);
            mRightMessages[currentLayer][i1] = Utils::BoxPlus(r1, l2 + r2);
            mRightMessages[currentLayer][i2] = Utils::BoxPlus(l1, r1) + r2;
        }
    }
}

bool Codec::PolarCode::BPDecoder::CheckStopCondition() const
{
    return mEncoder.Encode(GetInfVector()) == GetCodeword();
}

bool BPDecoder::CheckStopCondition(Codec::InfVector& codeword) const
{
    return mEncoder.Encode(GetInfVector()) == (codeword = GetCodeword());
}

template<typename T>
inline std::vector<T> BPDecoder::TransformLayer(
    size_t layer, std::function<T(double, double)> func
) const {
    std::vector<T> vec(mSpec->Length);
    for (size_t i = 0; i < mSpec->Length; i++) {
        vec[i] = func(mLeftMessages[layer][i], mRightMessages[layer][i]);
    }
    return vec;
}

std::vector<double> BPDecoder::GetLayerLLRs(size_t layer) const
{
    return TransformLayer<double>(layer, [](double l, double r) { return l + r; });
}

Codec::InfVector BPDecoder::GetLayerDecisions(size_t layer) const
{
    return TransformLayer<bool>(layer, [](double l, double r) { return l + r < 0; });
}

Codec::InfVector BPDecoder::GetInfVector() const
{
    return GetLayerDecisions(mSchedule[0]);
}

Codec::InfVector BPDecoder::GetCodeword() const
{
    return GetLayerDecisions(mNumLayers);
}
