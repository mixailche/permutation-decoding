#pragma once

#include "Construct/Perm.hpp"
#include "PolarSubcode/SCDecoder.h"
#include "Common.h"
#include "Specification.h"

namespace Codec {
    namespace PolarSubcode {
        class AutSCDecoder : public Decoder {
        public:
            using Perm = Construct::Perm;

            AutSCDecoder(
                const PolarSubcodeSpecification* spec,
                std::vector<Perm> permutations);

            ~AutSCDecoder() override = default;

            Codec::InfVector Decode(const std::vector<double>& inputLLRs) const;

        private:
            SCDecoder mDecoder;
            std::vector<Perm> mPermutations;
        };
    } // namespace PolarSubcode
} // namespace Codec
