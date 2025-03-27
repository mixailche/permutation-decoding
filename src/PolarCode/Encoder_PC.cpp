#include <algorithm>

#include "PolarCode/Encoder.h"

using Codec::PolarCodeSpecification;
using Codec::InfVector;

Codec::PolarCode::Encoder::Encoder(const PolarCodeSpecification* spec)
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

InfVector Codec::PolarCode::Encoder::Encode(const InfVector& infVector) const
{
    InfVector codeword(mSpec->Length, false);

    for (size_t i = 0, infPos = 0; i < mSpec->Length; i++) {
        if (mSpec->Frozen[i]) {
            codeword[i] = 0;
        }
        else {
            codeword[i] = infVector[infPos++];
        }
    }

    Polarize(codeword.begin(), codeword.end());

    return codeword;
}
