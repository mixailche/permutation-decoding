#pragma once

#include <vector>

#include "math/MatGF2.h"
#include "math/VecGF2.h"

namespace math {
    class Perm {
    public:
        Perm() = default;
        Perm(std::vector<size_t> newIndices);

        static Perm MakeIdentity(size_t length);
        static Perm MakeDigits(const std::vector<size_t>& digitsPerm);
        static Perm MakeAffine(const MatGF2& factor, const VecGF2& shift);

        size_t Length() const;

        size_t& operator[](size_t number);
        const size_t& operator[](size_t number) const;

        const std::vector<size_t>& AsVector() const;
        const MatGF2 AsMatGF2() const;

        template <typename T>
        std::vector<T> ApplyDirect(const std::vector<T>& vec) const;

        template <typename T>
        std::vector<T> ApplyReversed(const std::vector<T>& vec) const;

    private:
        std::vector<size_t> mNewIndices;
    };
} // namespace math

template <typename T>
std::vector<T> math::Perm::ApplyDirect(const std::vector<T>& vec) const
{
    std::vector<T> result(vec.size());

    for (size_t i = 0; i < vec.size(); i++) {
        result[this->operator[](i)] = vec[i];
    }

    return result;
}

template <typename T>
std::vector<T> math::Perm::ApplyReversed(const std::vector<T>& vec) const
{
    std::vector<T> result(vec.size());

    for (size_t i = 0; i < vec.size(); i++) {
        result[i] = vec[this->operator[](i)];
    }

    return result;
}
