#pragma once

#include <cstddef>
#include <vector>

#include "Specification.h"
#include "Common.h"

namespace Codec {
    namespace PolarCode {
        class Encoder : public Codec::Encoder {
        public:
            explicit Encoder(const PolarCodeSpecification* spec);
            ~Encoder() override = default;
            InfVector Encode(const InfVector& infVector) const override;

        private:
            const PolarCodeSpecification* mSpec;
        };
    } // namespace PolarCode
} // namespace Codec
