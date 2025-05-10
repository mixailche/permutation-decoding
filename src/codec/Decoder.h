#pragma once

#include <vector>

namespace codec {
    class Decoder {
    public:
        virtual std::vector<bool> Decode(const std::vector<double>& inputLLRs) const = 0;
        virtual size_t NumOperations() const = 0;
        virtual ~Decoder() = default;
    };
}