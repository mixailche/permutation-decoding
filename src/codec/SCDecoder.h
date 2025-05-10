#pragma once

#include "codec/Decoder.h"
#include "codec/PolarSpecification.h"

namespace codec {
    class SCDecoder : public Decoder {
    public:
        explicit SCDecoder(const PolarSpecification* spec);
        ~SCDecoder() override;

        std::vector<bool> Decode(const std::vector<double>& inputLLRs) const override;
        size_t NumOperations() const override;
        double Metric() const;

    private:
        const PolarSpecification* mSpec;
        size_t mNumLayers;

        mutable double mMetric;
        mutable size_t mNumOperations;
        mutable std::vector<double*> mLLRs;
        mutable std::vector<bool*> mSymbols;
        mutable std::vector<bool> mDynFrozenSymbols;

        inline size_t GetLayerSize(size_t layer) const;

        void ProcessLayer(size_t layer, size_t phase) const;
        void ProcessPhase(size_t phase) const;

        void AssignSymbol(size_t phase, bool value) const;
    };
}
