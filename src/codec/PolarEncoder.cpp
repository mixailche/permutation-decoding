#include <algorithm>
#include <vector>

#include "PolarEncoder.h"

using codec::PolarEncoder;
using codec::PolarSpecification;
using InfVector = std::vector<bool>;

PolarEncoder::PolarEncoder(const PolarSpecification* spec)
    : mSpec(spec)
{}

static void Polarize(InfVector::iterator begin, InfVector::iterator end)
{
    auto size = std::distance(begin, end);
    if (size == 1) {
        return;
    }

    auto half = size / 2;
    auto middle = begin + half;

    Polarize(begin, middle);
    Polarize(middle, end);

    std::transform(begin, middle, middle, begin, [](bool lhs, bool rhs) {
        return lhs ^ rhs;
    });
}

std::vector<bool> PolarEncoder::Encode(const std::vector<bool>& infVector) const
{
    std::vector<bool> codeword(mSpec->Length, false);

    for (size_t i = 0, infPos = 0; i < mSpec->Length; i++) {
        if (mSpec->StaticFrozen[i]) {
            codeword[i] = 0;
        }
        else if (!mSpec->Dynamic.Frozen[i]) {
            codeword[i] = infVector[infPos++];
        }

        for (auto dyn : mSpec->Dynamic.ForwardEquations[i]) {
            codeword[dyn] = codeword[dyn] ^ codeword[i];
        }
    }

    Polarize(codeword.begin(), codeword.end());

    return codeword;
}
