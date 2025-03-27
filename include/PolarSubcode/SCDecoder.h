#pragma once

#include "Common.h"

namespace Codec {
    namespace PolarSubcode {
        class SCDecoder : public Decoder {
        public:
            explicit SCDecoder(const PolarSubcodeSpecification* spec);
            ~SCDecoder() override;

            InfVector Decode(const std::vector<double>& inputLLRs) const override;
            double Metric() const;

        private:
            const PolarSubcodeSpecification* mSpec;
            size_t mNumLayers;

            mutable double mMetric;
            mutable std::vector<double*> mLLRs;
            mutable std::vector<bool*> mSymbols;
            mutable std::vector<bool> mDynFrozenSymbols;

            inline size_t GetLayerSize(size_t layer) const;

            void ProcessLayer(size_t layer, size_t phase) const;
            void ProcessPhase(size_t phase) const;

            void AssignSymbol(size_t phase, bool value) const;
        };
    } // namespace PolarSubcode
} // namespace Codec
