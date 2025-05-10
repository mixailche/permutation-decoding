#pragma once

#include <vector>
#include "PolarSpecification.h"

namespace codec {
    class PolarEncoder {
    public:
        explicit PolarEncoder(const PolarSpecification* spec);
        
        std::vector<bool> Encode(const std::vector<bool>& infVector) const;

    private:
        const PolarSpecification* mSpec;
    };
}
