#pragma once

#include <array>
#include <vector>
#include <cstddef>
#include <memory>
#include <any>

#include "Specification.h"

namespace Codec {
    using InfVector = std::vector<bool>;

    class Encoder {
    public:
        virtual ~Encoder() = default;

        virtual InfVector Encode(const InfVector& infVector) const = 0;
    };

    class Decoder {
    public:
        virtual ~Decoder() = default;
        virtual InfVector Decode(const std::vector<double>& inputLLRs) const = 0;
    };
} // namespace codec

namespace Running {    
    class ConsoleRunner {

    };

    class MatplotRunner {

    };
} // namespace running
