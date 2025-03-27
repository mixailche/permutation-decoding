#pragma once

#include "Construct/Perm.hpp"
#include "Common.h"
#include "Specification.h"

namespace Codec {
    namespace PolarCode {
        class AutBPDecoder : public Decoder {
        public:
            AutBPDecoder(
                const Codec::PolarCodeSpecification*, size_t numIterMax,
                std::vector<Perm> permutations
            );
            ~AutBPDecoder() override = default;

            InfVector Decode(const std::vector<double>& inputLLRs) const override;

        private:
            const Codec::PolarCodeSpecification* mSpec;
            std::vector<Perm> mPermutations;
            size_t mNumIterMax;
        };
    } // namespace PolarCode
} // namespace Codec
