#include "GField.h"

using math::GField;

GField::GField(uint64_t primitivePolynomialMask)
    : mNumDigits(0)
{
    size_t size = 1;

    for (auto copy = primitivePolynomialMask >> 1; copy > 0; copy >>= 1) {
        size <<= 1;
        mNumDigits++;
    }

    auto remainderMask = size - 1;
    auto reduction = primitivePolynomialMask & remainderMask;

    mLogarithms.resize(size);
    mExponents.resize(size - 1);
    Elem elem = 1;

    for (size_t deg = 0; deg < size - 1; deg++) {
        mLogarithms[elem] = deg;
        mExponents[deg] = elem;
        elem <<= 1;

        if ((elem & size) != 0) {
            elem = (elem & remainderMask) ^ reduction;
        }
    }
}

size_t GField::NumDigits() const noexcept
{
    return mNumDigits;
}

size_t GField::Size() const noexcept
{
    return mLogarithms.size();
}

size_t GField::Log(Elem elem) const
{
    return mLogarithms[elem];
}

GField::Elem GField::Exp(size_t deg) const
{
    return mExponents[deg];
}

GField::Elem GField::Pow(Elem lhs, Elem rhs) const
{
    if (lhs == 0) {
        return rhs == 0 ? 1 : 0;
    }

    auto currentDeg = mLogarithms[lhs];
    auto newDeg = (currentDeg * rhs) % (Size() - 1);
    return mExponents[newDeg];
}
