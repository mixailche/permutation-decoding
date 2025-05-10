#include <algorithm>

#include "math/MathUtils.h"
#include "codec/PermSCDecoder.h"
#include "codec/SCDecoder.h"

using codec::PermSCDecoder;

PermSCDecoder::PermSCDecoder(const PolarSpecification& spec, std::vector<Perm> permutations)
    : mPerms(std::move(permutations))
    , mPermutedSpecs(mPerms.size())
    , mNumOperations(0)
{
    auto freezingMatrix = spec.BuildFreezingMatrix();
    std::transform(
        mPerms.begin(), mPerms.end(), mPermutedSpecs.begin(),
        [&](const Perm& perm) {
            return math::TransformFreezingMatrix(freezingMatrix, perm.AsMatGF2());
        }
    );
}

size_t PermSCDecoder::NumOperations() const
{
    return mNumOperations;
}

std::vector<bool> PermSCDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    auto bestMetric = -std::numeric_limits<double>::max();
    std::vector<bool> bestCodeword;

    // metric evaluation -> n operations
    // metric comparison -> 1 operation
    mNumOperations = (inputLLRs.size() + 1) * mPerms.size();

    for (size_t i = 0; i < mPerms.size(); i++) {
        SCDecoder decoder(&mPermutedSpecs[i]);
        auto& perm = mPerms[i];
        auto permutedCodeword = decoder.Decode(perm.ApplyDirect(inputLLRs));
        auto codeword = perm.ApplyReversed(permutedCodeword);
        mNumOperations += decoder.NumOperations();

        double metric = 0;
        for (size_t i = 0; i < inputLLRs.size(); i++) {
            if ((inputLLRs[i] < 0) != codeword[i]) {
                metric -= std::abs(inputLLRs[i]);
            }
        }

        if (metric > bestMetric) {
            bestMetric = metric;
            bestCodeword = codeword;
        }
    }

    return bestCodeword;
}
