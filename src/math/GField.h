#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>

namespace math {
    class GField {
    public:
        GField(uint64_t primitivePolynomialMask);

        using Elem = uint64_t;

        size_t NumDigits() const noexcept;
        size_t Size() const noexcept;

        size_t Log(Elem elem) const;
        Elem Exp(size_t deg) const;
        Elem Pow(Elem lhs, Elem rhs) const;

    private:
        size_t mNumDigits;

        std::vector<size_t> mLogarithms;
        std::vector<Elem> mExponents;
    };
}
