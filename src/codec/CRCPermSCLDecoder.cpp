#include <algorithm>

#include "utils/Utils.hpp"
#include "math/MathUtils.h"
#include "codec/SCLDecoder.h"
#include "codec/CRCPermSCLDecoder.h"
#include <iostream>

using codec::CRCPermSCLDecoder;

CRCPermSCLDecoder::CRCPermSCLDecoder(const PolarSpecification& spec, std::vector<bool> generator, std::vector<Perm> perms, size_t maxPaths)
    : mGenerator(std::move(generator))
    , mPerms(std::move(perms))
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

std::vector<bool> CRCPermSCLDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    auto length = inputLLRs.size();
    auto dimension = mPermutedSpecs[0].Dimension;
    auto kernel = math::BuildArikanKernel(utils::IntLog2(length)).Transpose();
    auto bestMetric = -std::numeric_limits<double>::max();
    std::vector<bool> bestCodeword;

    // metric evaluation -> n operations
    // metric comparison -> 1 operation
    mNumOperations = (inputLLRs.size() + 1) * mPerms.size();
    auto numCRCBits = mGenerator.size() - 1;
    
    for (size_t i = 0; i < mPerms.size(); i++) {
        SCLDecoder decoder(&mPermutedSpecs[i], mMaxPaths);
        auto& perm = mPerms[i];
        auto permutedCodeword = decoder.Decode(perm.ApplyDirect(inputLLRs));
        auto codeword = perm.ApplyReversed(permutedCodeword);
        auto u = kernel * codeword;
        
        std::vector<bool> infVector(dimension);
        for (size_t i = 0, j = 0; i < length; i++) {
            if (!mPermutedSpecs[0].StaticFrozen[i]) {
                infVector[j++] = u[i];
            }
        }

        auto middle = infVector.begin() + dimension - numCRCBits;
        std::vector<bool> payload(infVector.begin(), middle);
        std::vector<bool> crc(middle, infVector.end());

        if (math::CalculateCRC(payload, mGenerator) != crc) {
            continue;
        }

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

size_t CRCPermSCLDecoder::NumOperations() const
{
    return mNumOperations;
}
