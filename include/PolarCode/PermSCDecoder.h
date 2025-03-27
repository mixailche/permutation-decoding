#pragma once

#include "PolarCode/SCDecoder.h"
#include "Construct/Perm.hpp"
#include "Specification.h"
#include "Common.h"

namespace Codec {
    namespace PolarCode {
        class PermSCDecoder : public Decoder {
        public:
            using Perm = Construct::Perm;

            PermSCDecoder(const PolarCodeSpecification* spec, std::vector<Perm> permutations);

            ~PermSCDecoder() override = default;

            Codec::InfVector Decode(const std::vector<double>& inputLLRs) const;

        private:
            const std::vector<bool>& mFrozen;
            std::vector<Construct::Perm> mPerms;
        };
    } // namespace PolarCode
} // namespace Codec
