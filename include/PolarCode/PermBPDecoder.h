#pragma once

#include "Construct/Perm.hpp"
#include "Common.h"

namespace Codec {
    namespace PolarCode {
        class PermBPDecoder : public Decoder {
        public:
            PermBPDecoder(
                const Codec::PolarCodeSpecification* spec, size_t numIterMax,
                std::vector<Perm> permutaions
            );
            ~PermBPDecoder() override = default;

            InfVector Decode(const std::vector<double>& inputLLRs) const override;

        private:
            const Codec::PolarCodeSpecification* mSpec;
            std::vector<Perm> mPermutations;
            size_t mNumIterMax;
        };
    } // namespace PolarCode
} // namespace Codec

