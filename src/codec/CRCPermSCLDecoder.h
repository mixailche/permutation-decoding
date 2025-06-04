#pragma once

#include <vector>

#include "math/Perm.hpp"
#include "codec/PolarSpecification.h"
#include "codec/Decoder.h"

namespace codec {
    class CRCPermSCLDecoder : public Decoder {
    public:
        using Perm = math::Perm;

        CRCPermSCLDecoder(
            const PolarSpecification& spec,
            std::vector<bool> generator, std::vector<Perm> perms, size_t maxPaths);
        
        ~CRCPermSCLDecoder() override = default;

        std::vector<bool> Decode(const std::vector<double>& inputLLRs) const override;
        size_t NumOperations() const override;

    private:
        std::vector<bool> mGenerator;
        std::vector<Perm> mPerms;
        std::vector<PolarSpecification> mPermutedSpecs;
        size_t mMaxPaths;
        mutable size_t mNumOperations;
    };
}
