#include <algorithm>

#include "math/MathUtils.h"
#include "codec/SCLDecoder.h"
#include "codec/PermSCLDecoder.h"
#include <iostream>

using codec::PermSCLDecoder;

PermSCLDecoder::PermSCLDecoder(const PolarSpecification& spec, std::vector<Perm> perms, size_t maxPaths)
    : mPerms(std::move(perms))
    , mPermutedSpecs(mPerms.size())
    , mMaxPaths(maxPaths)
    , mNumOperations(0)
{
    auto freezingMatrix = spec.BuildFreezingMatrix();
    std::transform(
        mPerms.begin(), mPerms.end(), mPermutedSpecs.begin(),
        [&](const Perm& perm) -> PolarSpecification {
            return math::TransformFreezingMatrix(freezingMatrix, perm.AsMatGF2());
        }
    );
}

std::vector<bool> PermSCLDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    auto bestMetric = -std::numeric_limits<double>::max();
    std::vector<bool> bestCodeword;

    // metric evaluation -> n operations
    // metric comparison -> 1 operation
    mNumOperations = (inputLLRs.size() + 1) * mPerms.size();

    for (size_t i = 0; i < mPerms.size(); i++) {
        SCLDecoder decoder(&mPermutedSpecs[i], mMaxPaths);
        auto& perm = mPerms[i];
        auto permutedCodeword = decoder.Decode(perm.ApplyReversed(inputLLRs));
        auto codeword = perm.ApplyDirect(permutedCodeword);
        mNumOperations += decoder.NumOperations();

        double metric = 0;
        for (size_t j = 0; j < inputLLRs.size(); j++) {
            if ((inputLLRs[j] < 0) != codeword[j]) {
                metric -= std::abs(inputLLRs[j]);
            }
        }

        if (metric > bestMetric) {
            bestMetric = metric;
            bestCodeword = codeword;
        }
    }

    return bestCodeword;
}

size_t PermSCLDecoder::NumOperations() const
{
    return mNumOperations;
}
