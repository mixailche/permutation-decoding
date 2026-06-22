#pragma once

#include <vector>

#include "codec/Decoder.h"
#include "codec/PolarSpecification.h"

namespace codec {
    class ConstituentPermSCLDecoder : public Decoder {
    public:
        explicit ConstituentPermSCLDecoder(const PolarSpecification* spec, size_t numPerms, size_t minNumECs);
        ~ConstituentPermSCLDecoder() override;

        std::vector<bool> Decode(const std::vector<double>& inputLLRs) const override;
        size_t NumOperations() const override;

    private:
        const PolarSpecification* mSpec;
        size_t mNumLayers;
        std::vector<std::vector<Decoder*>> mOuterDecoders;
        std::vector<PolarSpecification*> mOuterSpecs;

        mutable size_t mNumOperations;
        mutable std::vector<double*> mLLRs;
        mutable std::vector<bool*> mSymbols;
        mutable std::vector<bool> mDynFrozenSymbols;

        void InitOuterDecoders(size_t numPerms, size_t minNumECs, size_t layer, size_t phase);

        inline size_t GetLayerSize(size_t layer) const;

        void ProcessNode(size_t layer, size_t phase) const;
        void ApplyOuterDecoder(size_t layer, size_t phase) const;
    };
};
