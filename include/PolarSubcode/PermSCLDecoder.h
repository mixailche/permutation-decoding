#pragma once

#include "Common.h"
#include "Construct/Perm.hpp"

namespace Codec {
    namespace PolarSubcode {
        class PermSCLDecoder : public Decoder {
        public:
            using Perm = Construct::Perm;

            PermSCLDecoder(
                const PolarSubcodeSpecification* spec,
                std::vector<Perm> permutations,
                size_t maxPaths
            );

            ~PermSCLDecoder() override = default;

            Codec::InfVector Decode(const std::vector<double>& inputLLRs) const;

        private:
            std::vector<Perm> mPermutations;
            std::vector<Codec::PolarSubcodeSpecification> mPermutedSpecs;

            size_t mMaxPaths;
        };
    } // namespace PolarCode
} // namespace Codec

