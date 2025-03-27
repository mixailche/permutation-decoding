#include <numeric>

#include "PolarSubcode/AutSCDecoer.h"
#include "PolarSubcode/SCDecoder.h"

using Codec::PolarSubcode::AutSCDecoder;

AutSCDecoder::AutSCDecoder(
    const PolarSubcodeSpecification* spec,
    std::vector<Perm> permutations)
    : mDecoder(spec)
    , mPermutations(std::move(permutations))
{
}

Codec::InfVector AutSCDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    auto bestMetric = -std::numeric_limits<double>::max();
    Codec::InfVector bestCodeword;

    for (auto& perm : mPermutations) {
        auto permutedCodeword = mDecoder.Decode(perm.ApplyDirect(inputLLRs));
        auto codeword = perm.ApplyReversed(permutedCodeword);
        
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
