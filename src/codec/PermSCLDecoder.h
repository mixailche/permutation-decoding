#pragma once

#include <vector>

#include "math/Perm.hpp"
#include "codec/PolarSpecification.h"
#include "codec/Decoder.h"

namespace codec {
    class PermSCLDecoder : public Decoder {
    public:
        using Perm = math::Perm;

        PermSCLDecoder(const PolarSpecification& spec, std::vector<Perm> perms, size_t maxPaths);
        ~PermSCLDecoder() override = default;

        std::vector<bool> Decode(const std::vector<double>& inputLLRs) const;
        size_t NumOperations() const;

    private:
        std::vector<Perm> mPerms;
        std::vector<PolarSpecification> mPermutedSpecs;
        size_t mMaxPaths;
        mutable size_t mNumOperations;
    };
}
