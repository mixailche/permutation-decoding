#include "PolarCode/BPDecoder.h"
#include "PolarCode/AutBPDecoder.h"

using Codec::PolarCode::AutBPDecoder;

AutBPDecoder::AutBPDecoder(
    const Codec::PolarCodeSpecification* spec, size_t numIterMax,
    std::vector<Perm> permutations)
    : mSpec(spec)
    , mPermutations(std::move(permutations))
    , mNumIterMax(numIterMax)
{}

Codec::InfVector AutBPDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    auto bestMetric = std::numeric_limits<double>::max();
    Codec::InfVector bestCodeword;
 
    for (const auto& perm : mPermutations) {
        BPDecoder decoder(mSpec, mNumIterMax);
        auto codeword = perm.ApplyReversed(decoder.Decode(perm.ApplyDirect(inputLLRs)));

        double metric = 0;
        for (size_t i = 0; i < codeword.size(); i++) {
            if (codeword[i] != (inputLLRs[i] < 0)) {
                metric += std::abs(inputLLRs[i]);
            }
        }

        if (metric < bestMetric) {
            bestMetric = metric;
            bestCodeword = std::move(codeword);
        }
    }

    return bestCodeword;
}
