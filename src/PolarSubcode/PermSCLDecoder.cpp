#include <algorithm>

#include "PolarSubcode/SCLDecoder.h"
#include "PolarSubcode/PermSCLDecoder.h"

using Codec::PolarSubcode::PermSCLDecoder;

PermSCLDecoder::PermSCLDecoder(
    const PolarSubcodeSpecification* spec,
    std::vector<Perm> permutations, size_t maxPaths)
    : mPermutations(std::move(permutations))
    , mPermutedSpecs(mPermutations.size())
    , mMaxPaths(maxPaths)
{
    std::transform(
        mPermutations.begin(), mPermutations.end(), mPermutedSpecs.begin(),
        [&](const Perm& perm) {
            return Construct::PermuteSpecification(*spec, perm);
        }
    );
}

Codec::InfVector PermSCLDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    auto bestMetric = -std::numeric_limits<double>::max();
    Codec::InfVector bestCodeword;
    bool diff = false;

    for (size_t i = 0; i < mPermutations.size(); i++) {
        SCLDecoder scl(&mPermutedSpecs[i], mMaxPaths);
        auto& perm = mPermutations[i];
        auto permutedCodeword = scl.Decode(perm.ApplyReversed(inputLLRs));
        auto codeword = perm.ApplyDirect(permutedCodeword);

        double metric = 0;
        for (size_t j = 0; j < inputLLRs.size(); j++) {
            if ((inputLLRs[j] < 0) != codeword[j]) {
                metric -= std::abs(inputLLRs[j]);
            }
        }

        if (metric > bestMetric) {
            if (i != 0 && !diff) {
                diff = true;
            }
            bestMetric = metric;
            bestCodeword = codeword;
        }
    }

    return bestCodeword;
}