#pragma once

#include "Common.h"
#include "Specification.h"

namespace Codec {
    namespace PolarSubcode {
        class Encoder : public Codec::Encoder {
        public:
            explicit Encoder(const PolarSubcodeSpecification* spec);
            ~Encoder() override = default;
            InfVector Encode(const InfVector& infVector) const override;

        private:
            const PolarSubcodeSpecification* mSpec;
        };
    } // namespace PolarSubcode
} // namespace Codec
