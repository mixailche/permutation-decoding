#include "PolarCode/BPPermSCDecoder.h"
#include "PolarCode/BPDecoder.h"
#include "PolarCode/SCDecoder.h"

using Codec::PolarCode::BPPermSCDecoder;

BPPermSCDecoder::BPPermSCDecoder(
    const PolarCodeSpecification* spec, size_t numIterMax,
    std::vector<Perm> permutations)
    : mSpec(spec)
    , mPerms(std::move(permutations))
    , mNumIterMax(numIterMax)
{
}

Codec::InfVector BPPermSCDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    BPDecoder bp(mSpec, mNumIterMax);
    auto [success, llrs] = bp.SoftDecode(inputLLRs);
    
    if (success) {
        std::vector<bool> codeword(mSpec->Length);
        std::transform(llrs.begin(), llrs.end(), codeword.begin(),
            [](double y) { return y < 0; });
        return codeword;
    }

    auto bestMetric = -std::numeric_limits<double>::max();
    Codec::InfVector bestCodeword;

    for (const auto& perm : mPerms) {
        PolarCodeSpecification spec = perm.ApplyDirect(mSpec->Frozen);
        SCDecoder sc(&spec);

        auto permutedCodeword = sc.Decode(perm.ApplyDirect(llrs));
        auto metric = sc.Metric();

        if (metric > bestMetric) {
            bestMetric = metric;
            bestCodeword = perm.ApplyReversed(permutedCodeword);
        }
    }

    return bestCodeword;
}

