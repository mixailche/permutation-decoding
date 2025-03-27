#pragma once

#include <vector>
#include "NTL/mat_GF2.h"

namespace Codec {
    struct PolarCodeSpecification {
        std::vector<bool> Frozen;
        size_t Length;
        size_t Dimension;

        PolarCodeSpecification(std::vector<bool> frozen);
    };

    struct PolarSubcodeSpecification {
        PolarSubcodeSpecification();
        PolarSubcodeSpecification(const PolarSubcodeSpecification& other);
        PolarSubcodeSpecification(PolarSubcodeSpecification&& other);

        PolarSubcodeSpecification& operator=(const PolarSubcodeSpecification& other);
        PolarSubcodeSpecification& operator=(PolarSubcodeSpecification&& other);

        explicit PolarSubcodeSpecification(size_t length);
        PolarSubcodeSpecification(const NTL::mat_GF2& freezingMatrix);

        std::vector<bool> StaticFrozen;

        struct {
            std::vector<bool> Frozen;
            std::vector<std::vector<size_t>> ForwardEquations;
        } Dynamic;

        size_t Length;
        size_t Dimension;
    private:
        PolarSubcodeSpecification(size_t length, size_t dimension);
    };
}
