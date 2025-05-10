#include <algorithm>

#include "math/MathUtils.h"
#include "codec/PolarSpecification.h"

using codec::PolarSpecification;

PolarSpecification::PolarSpecification(size_t length, size_t dimension)
    : StaticFrozen(length, false)
    , Dynamic({
        std::vector(length, false),
        std::vector(length, std::vector<size_t>())
    })
    , Length(length)
    , Dimension(dimension)
{}

PolarSpecification::PolarSpecification()
    : PolarSpecification(0, 0)
{}

PolarSpecification::PolarSpecification(const PolarSpecification& other)
    : StaticFrozen(other.StaticFrozen)
    , Dynamic(other.Dynamic)
    , Length(other.Length)
    , Dimension(other.Dimension)
{}

PolarSpecification::PolarSpecification(PolarSpecification&& other) noexcept
    : StaticFrozen(std::move(other.StaticFrozen))
    , Dynamic({
        std::move(other.Dynamic.Frozen),
        std::move(other.Dynamic.ForwardEquations)
        })
    , Length(other.Length)
    , Dimension(other.Dimension)
{}

PolarSpecification& PolarSpecification::operator=(const PolarSpecification& other)
{
    StaticFrozen = other.StaticFrozen;
    Dynamic = other.Dynamic;
    Length = other.Length;
    Dimension = other.Dimension;
    return *this;
}

PolarSpecification& PolarSpecification::operator=(PolarSpecification&& other) noexcept
{
    StaticFrozen = std::move(other.StaticFrozen);
    Dynamic.Frozen = std::move(other.Dynamic.Frozen);
    Dynamic.ForwardEquations = std::move(other.Dynamic.ForwardEquations);
    Length = other.Length;
    Dimension = other.Dimension;
    return *this;
}

PolarSpecification::PolarSpecification(size_t length)
    : PolarSpecification(length, length)
{}

PolarSpecification::PolarSpecification(const math::MatGF2& freezingMatrix)
    : PolarSpecification(freezingMatrix.NumCols())
{
    std::vector<std::vector<size_t>> equations(Length);

    for (size_t r = 0; r < freezingMatrix.NumRows(); r++) {
        auto& row = freezingMatrix[r];
        if (auto end = row.End()) {
            auto i = *end;
            auto& equation = equations[i];

            for (size_t j = 0; j < i; j++) {
                if (row[j]) {
                    equation.push_back(j);
                }
            }

            if (equation.empty()) {
                StaticFrozen[i] = true;
            }
        }
    }

    for (size_t i = 0; i < Length; i++) {
        for (size_t j : equations[i]) {
            if (!StaticFrozen[j]) {
                Dynamic.Frozen[i] = true;
                Dynamic.ForwardEquations[j].push_back(i);
            }
        }

        if (!equations[i].empty() && !Dynamic.Frozen[i]) {
            // all summands are trivial
            StaticFrozen[i] = true;
        }
    }

    Dimension =
        std::count(StaticFrozen.begin(), StaticFrozen.end(), 0) -
        std::count(Dynamic.Frozen.begin(), Dynamic.Frozen.end(), 1);
}

math::MatGF2 PolarSpecification::BuildFreezingMatrix() const
{
    auto freezingMatrix = math::MatGF2::Zeros(Length - Dimension, Length);
    std::vector<size_t> dynamicFrozenSymbolRows(Length);

    for (size_t i = 0, j = 0; j < Length; i++, j++) {
        if (StaticFrozen[j]) {
            freezingMatrix.Set(i, j, 1);
        }
        else if (Dynamic.Frozen[j]) {
            freezingMatrix.Set(i, j, 1);
            dynamicFrozenSymbolRows[j] = i;
        }
        else {
            i--;
        }
    }

    for (size_t i = 0; i < Length; i++) {
        for (size_t dyn : Dynamic.ForwardEquations[i]) {
            freezingMatrix.Set(dynamicFrozenSymbolRows[dyn], i, 1);
        }
    }

    return freezingMatrix;
}
