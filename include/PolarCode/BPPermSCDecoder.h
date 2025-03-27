#pragma once

#include "Common.h"
#include "Specification.h"
#include "Construct/Perm.hpp"

namespace Codec {
    namespace PolarCode {
        class BPPermSCDecoder : public Decoder {
        public:
            using Perm = Construct::Perm;

            BPPermSCDecoder(
                const PolarCodeSpecification* spec, size_t numIterMax,
                std::vector<Perm> permutations
            );

            ~BPPermSCDecoder() override = default;

            Codec::InfVector Decode(const std::vector<double>& inputLLRs) const;

        private:
            const PolarCodeSpecification* mSpec;
            std::vector<Construct::Perm> mPerms;
            size_t mNumIterMax;
        };
    } // namespace PolarCode
} // namespace Codec

