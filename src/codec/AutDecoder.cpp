#include "codec/AutDecoder.h"

using codec::AutDecoder;

AutDecoder::AutDecoder(const Decoder* decoder, std::vector<math::Perm> perms)
    : mDecoder(decoder)
    , mNumOperations(0)
    , mPerms(std::move(perms))
{}

std::vector<bool> AutDecoder::Decode(const std::vector<double>& inputLLRs) const
{
    auto bestMetric = -std::numeric_limits<double>::max();
    std::vector<bool> bestCodeword;

    // metric evaluation -> n operations
    // metric comparison -> 1 operation
    mNumOperations = (inputLLRs.size() + 1) * mPerms.size();

    for (auto& perm : mPerms) {
        auto permutedCodeword = mDecoder->Decode(perm.ApplyDirect(inputLLRs));
        auto codeword = perm.ApplyReversed(permutedCodeword);
        mNumOperations += mDecoder->NumOperations();

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

size_t codec::AutDecoder::NumOperations() const
{
    return mNumOperations;
}
