#pragma once

#include <memory>

#include "PolarSubcode/SCDecoder.h"
#include "Construct/Perm.hpp"
#include "Specification.h"
#include "Common.h"

namespace Codec {
    namespace PolarSubcode {
        class PermSCDecoder : public Decoder {
        public:
            using Perm = Construct::Perm;

            PermSCDecoder(
                const PolarSubcodeSpecification* spec, 
                std::vector<Perm> permutations
            );

            ~PermSCDecoder() override = default;

            Codec::InfVector Decode(const std::vector<double>& inputLLRs) const;

        private:
            std::vector<Perm> mPermutations;
            std::vector<Codec::PolarSubcodeSpecification> mPermutedSpecs;
        };
    } // namespace PolarCode
} // namespace Codec
