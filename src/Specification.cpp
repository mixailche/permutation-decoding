#include <algorithm>
#include "NTLHelp.h"
#include "Specification.h"

using Codec::PolarCodeSpecification;
using Codec::PolarSubcodeSpecification;

PolarCodeSpecification::PolarCodeSpecification(std::vector<bool> frozen)
    : Frozen(std::move(frozen))
    , Length(Frozen.size())
    , Dimension(std::count(Frozen.begin(), Frozen.end(), false))
{}

PolarSubcodeSpecification::PolarSubcodeSpecification(size_t length, size_t dimension)
    : StaticFrozen(length, false)
    , Dynamic({
        std::vector(length, false),
        std::vector(length, std::vector<size_t>())
    })
    , Length(length)
    , Dimension(dimension)
{}

PolarSubcodeSpecification::PolarSubcodeSpecification()
    : PolarSubcodeSpecification(0, 0)
{}

PolarSubcodeSpecification::PolarSubcodeSpecification(const PolarSubcodeSpecification& other)
    : StaticFrozen(other.StaticFrozen)
    , Dynamic(other.Dynamic)
    , Length(other.Length)
    , Dimension(other.Dimension)
{}

PolarSubcodeSpecification::PolarSubcodeSpecification(PolarSubcodeSpecification&& other)
    : StaticFrozen(std::move(other.StaticFrozen))
    , Dynamic({
        std::move(other.Dynamic.Frozen),
        std::move(other.Dynamic.ForwardEquations)
    })
    , Length(other.Length)
    , Dimension(other.Dimension)
{}

PolarSubcodeSpecification& Codec::PolarSubcodeSpecification::operator=(const PolarSubcodeSpecification & other)
{
    StaticFrozen = other.StaticFrozen;
    Dynamic = other.Dynamic;
    Length = other.Length;
    Dimension = other.Dimension;
    return *this;
}

PolarSubcodeSpecification& Codec::PolarSubcodeSpecification::operator=(PolarSubcodeSpecification&& other)
{
    StaticFrozen = std::move(other.StaticFrozen);
    Dynamic.Frozen = std::move(other.Dynamic.Frozen);
    Dynamic.ForwardEquations = std::move(other.Dynamic.ForwardEquations);
    Length = other.Length;
    Dimension = other.Dimension;
    return *this;
}

PolarSubcodeSpecification::PolarSubcodeSpecification(size_t length)
    : PolarSubcodeSpecification(length, length)
{}

PolarSubcodeSpecification::PolarSubcodeSpecification(const NTL::mat_GF2& freezingMatrix)
    : PolarSubcodeSpecification(freezingMatrix.NumCols())
{
    std::vector<std::vector<size_t>> equations(Length);

    for (size_t r = 0; r < freezingMatrix.NumRows(); r++) {
        auto& row = freezingMatrix[r];
        if (auto end = Math::HighestOnePos(row)) {
            auto i = *end;
            auto& equation = equations[i];

            for (size_t j = 0; j < i; j++) {
                if (NTL::IsOne(row[j])) {
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
