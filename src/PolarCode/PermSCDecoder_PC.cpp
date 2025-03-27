#include "PolarCode/PermSCDecoder.h"
#include "PolarCode/SCDecoder.h"

using Codec::PolarCode::PermSCDecoder;
using Codec::PolarCode::SCDecoder;

PermSCDecoder::PermSCDecoder(const PolarCodeSpecification* spec, std::vector<Perm> permutations)
    : mFrozen(spec->Frozen)
    , mPerms(std::move(permutations))
{}

Codec::InfVector PermSCDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    auto bestMetric = -std::numeric_limits<double>::max();
    Codec::InfVector bestCodeword;
    size_t wrong = 0;
    for (const auto& perm : mPerms) {
        /*if (mFrozen != perm.ApplyDirect(mFrozen)) {
            std::cout << wrong++ << std::endl;
        }*/
        //PolarCodeSpecification spec = perm.ApplyDirect(mFrozen);
        PolarCodeSpecification spec = mFrozen;
        SCDecoder sc(&spec);
        
        auto permutedCodeword = sc.Decode(perm.ApplyDirect(inputLLRs));
        auto metric = sc.Metric();

        if (metric > bestMetric) {
            bestMetric = metric;
            bestCodeword = perm.ApplyReversed(permutedCodeword);
        }
    }

    return bestCodeword;
}
