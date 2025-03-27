#pragma once

#include <cstddef>
#include <vector>

#include "Utils.hpp"
#include "Common.h"
#include "Specification.h"

namespace Codec {
    namespace PolarCode {
        class SCDecoder : public Decoder {
        public:
            explicit SCDecoder(const PolarCodeSpecification* spec);
            ~SCDecoder() override;

            InfVector Decode(const std::vector<double>& inputLLRs) const override;
            double Metric() const;

        private:
            const PolarCodeSpecification* mSpec;
            size_t mNLayers;

            mutable double mMetric;
            mutable std::vector<double*> mLLRs;
            mutable std::vector<bool*> mSymbols;
        
            inline size_t GetLayerSize(size_t layer) const;

            void ProcessLayer(size_t layer, size_t phase) const;
            void ProcessPhase(size_t phase) const;
        };
    } // namespace PolarCode
} // namespace Codec
