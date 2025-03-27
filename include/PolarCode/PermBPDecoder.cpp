#include <algorithm>

#include "Construct/Perm.hpp"
#include "PolarCode/BPDecoder.h"
#include "PolarCode/PermBPDecoder.h"

using Codec::PolarCode::PermBPDecoder;

PermBPDecoder::PermBPDecoder(
    const Codec::PolarCodeSpecification* spec, size_t numIterMax,
    std::vector<Perm> permutations)
    : mSpec(spec)
    , mPermutations(std::move(permutations))
    , mNumIterMax(numIterMax)
{
}

Codec::InfVector PermBPDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    //auto bestMetric = std::numeric_limits<double>::max();
    //Codec::InfVector bestCodeword;

    std::vector<double> llrs;

    for (const auto& perm : mPermutations) {
        Codec::PolarCodeSpecification permutedSpec = perm.ApplyDirect(mSpec->Frozen);
        BPDecoder bp(&permutedSpec, mNumIterMax);

        auto [success, decoded] = bp.SoftDecode(perm.ApplyDirect(
            llrs.empty() ? inputLLRs : llrs
        ));
        llrs = std::move(perm.ApplyReversed(decoded));
        
        if (success) {
            break;
        }

        /*auto decoded = bp.Decode(perm.ApplyDirect(llrs));
        auto codeword = perm.ApplyReversed(decoded);

        double metric = 0;
        for (size_t i = 0; i < codeword.size(); i++) {
            if (codeword[i] != (inputLLRs[i] < 0)) {
                metric += std::abs(inputLLRs[i]);
            }
        }

        if (metric < bestMetric) {
            bestMetric = metric;
            bestCodeword = std::move(codeword);
        }*/
    }
    
    std::vector<bool> codeword(mSpec->Length);
    std::transform(llrs.begin(), llrs.end(), codeword.begin(), [](double y) { return y < 0; });
    
    return codeword;
}


