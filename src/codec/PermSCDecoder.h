#pragma once

#include <vector>

#include "math/Perm.hpp"
#include "codec/PolarSpecification.h"
#include "codec/Decoder.h"

namespace codec {
    class PermSCDecoder : public Decoder {
    public:
        using Perm = math::Perm;

        PermSCDecoder(const PolarSpecification& spec, std::vector<Perm> permutations);
        ~PermSCDecoder() override = default;

        std::vector<bool> Decode(const std::vector<double>&inputLLRs) const;
        size_t NumOperations() const;

    private:
        std::vector<Perm> mPerms;
        std::vector<PolarSpecification> mPermutedSpecs;
        mutable size_t mNumOperations;
    };
}
