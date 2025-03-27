#include "Linalg/Field.h"

using Linalg::Field;

static std::vector<std::vector<size_t>> BuildCyclotomicCosets(size_t size)
{
    std::vector used(size - 1, false);
    std::vector<std::vector<size_t>> result = { { 0 } };

    for (size_t leader = 1; leader < size - 1; leader++) {
        if (used[leader]) {
            continue;
        }
        std::vector<size_t> coset = { leader };
        
        for (size_t deg = leader * 2; deg != leader; deg = (deg * 2) % (size - 1)) {
            coset.push_back(deg);
        }

        for (auto element : coset) {
            used[element] = true;
        }

        result.push_back(coset);
    }

    return result;
}

Field::Field(uint64_t primitivePolynomialMask)
    : mNDigits(0)
{
    size_t size = 1;

    for (auto copy = primitivePolynomialMask >> 1; copy > 0; copy >>= 1) {
        size <<= 1;
        mNDigits++;
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

    mCyclotomicCosets = BuildCyclotomicCosets(size);
    mCosetLeaders.resize(size - 1);

    for (const auto& coset : mCyclotomicCosets) {
        auto leader = coset.front();
        for (auto element : coset) {
            mCosetLeaders[element] = leader;
        }
    }
}

size_t Field::NDigits() const noexcept
{
    return mNDigits;
}

size_t Field::Size() const noexcept
{
    return mLogarithms.size();
}

size_t Field::Log(Elem elem) const
{
    return mLogarithms[elem];
}

Field::Elem Field::Exp(size_t deg) const
{
    return mExponents[deg];
}

Field::Elem Field::Pow(Elem lhs, Elem rhs) const
{
    if (lhs == 0) {
        return rhs == 0 ? 1 : 0;
    }

    auto currentDeg = mLogarithms[lhs];
    auto newDeg = (currentDeg * rhs) % (Size() - 1);
    return mExponents[newDeg];
}

const std::vector<std::vector<size_t>>& Field::CyclotomicCosets() const
{
    return mCyclotomicCosets;
}

const std::vector<size_t>& Field::CyclotomicCoset(size_t deg) const
{
    return mCyclotomicCosets[mCosetLeaders[deg]];
}


