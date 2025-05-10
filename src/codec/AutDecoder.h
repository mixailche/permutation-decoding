#pragma once

#include "math/Perm.hpp"
#include "codec/Decoder.h"

namespace codec {
    class AutDecoder : public Decoder {
    public:
        AutDecoder(const Decoder* decoder, std::vector<math::Perm> perms);
        ~AutDecoder() override = default;

        std::vector<bool> Decode(const std::vector<double>& inputLLRs) const override;
        size_t NumOperations() const override;

    private:
        const Decoder* mDecoder;
        mutable size_t mNumOperations;
        std::vector<math::Perm> mPerms;
    };
}
