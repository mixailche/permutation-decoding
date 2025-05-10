#pragma once

#include "math/Perm.hpp"
#include "math/MatGF2.h"

namespace codec {
    class PolarSpecification {
    public:
        PolarSpecification();
        PolarSpecification(const PolarSpecification& other);
        PolarSpecification(PolarSpecification&& other) noexcept;

        PolarSpecification& operator=(const PolarSpecification& other) ;
        PolarSpecification& operator=(PolarSpecification&& other) noexcept;

        explicit PolarSpecification(size_t length);
        PolarSpecification(const math::MatGF2& freezingMatrix);

        std::vector<bool> StaticFrozen;

        struct {
            std::vector<bool> Frozen;
            std::vector<std::vector<size_t>> ForwardEquations;
        } Dynamic;

        size_t Length;
        size_t Dimension;

        math::MatGF2 BuildFreezingMatrix() const;

    private:
        PolarSpecification(size_t length, size_t dimension);
    };
}