#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>

namespace Linalg {
    class Field {
    public:
        Field(uint64_t primitivePolynomialMask);

        using Elem = uint64_t;

        size_t NDigits() const noexcept;
        size_t Size() const noexcept;

        size_t Log(Elem elem) const;
        Elem Exp(size_t deg) const;
        Elem Pow(Elem lhs, Elem rhs) const;

        const std::vector<std::vector<size_t>>& CyclotomicCosets() const;
        const std::vector<size_t>& CyclotomicCoset(size_t deg) const;

    private:
        size_t mNDigits;
        
        std::vector<size_t> mLogarithms;
        std::vector<Elem> mExponents;

        std::vector<std::vector<size_t>> mCyclotomicCosets;
        std::vector<size_t> mCosetLeaders;
    };
}
