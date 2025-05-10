#include "utils/Utils.hpp"
#include "math/MathUtils.h"
#include "math/Perm.hpp"

using math::Perm;
using math::MatGF2;
using math::VecGF2;

Perm::Perm(std::vector<size_t> newIndices)
    : mNewIndices(std::move(newIndices))
{}

Perm Perm::MakeIdentity(size_t length)
{
    return utils::Iota(length);
}

static size_t ApplyDigitsPerm(const std::vector<size_t>& perm, size_t number)
{
    size_t result = 0;
    for (size_t i = 0; i < perm.size(); i++) {
        if (number & (1ull << perm[i])) {
            result |= (1ull << i);
        }
    }
    return result;
}

Perm Perm::MakeDigits(const std::vector<size_t>& digitsPerm)
{
    std::vector<size_t> newIndices(1ull << digitsPerm.size());
    for (size_t i = 0; i < newIndices.size(); i++) {
        newIndices[i] = ApplyDigitsPerm(digitsPerm, i);
    }
    return std::move(newIndices);
}

Perm Perm::MakeAffine(const MatGF2& factor, const VecGF2& shift)
{
    size_t length = factor.NumRows();
    size_t cardinality = 1ull << length;
    std::vector<size_t> newIndices(cardinality);

    for (size_t i = 0; i < cardinality; i++) {
        auto digitsVector = math::IndexToVecGF2(i, length);
        newIndices[i] = math::VecGF2ToIndex(factor * digitsVector + shift);
    }

    return newIndices;
}

size_t Perm::Length() const
{
    return mNewIndices.size();
}

size_t& Perm::operator[](size_t number)
{
    return mNewIndices[number];
}

const size_t& Perm::operator[](size_t number) const
{
    return mNewIndices[number];
}

const std::vector<size_t>& Perm::AsVector() const
{
    return mNewIndices;
}

const MatGF2 Perm::AsMatGF2() const
{
    auto permMatrix = MatGF2::Zeros(Length(), Length());
    for (size_t i = 0; i < mNewIndices.size(); i++) {
        permMatrix.Set(i, mNewIndices[i], 1);
    }
    return permMatrix;
}
