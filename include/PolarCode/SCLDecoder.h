#pragma once

#include <stack>

#include "Specification.h"
#include "Common.h"

namespace Codec {
    namespace PolarCode {
        class SCLDecoder : public Decoder {
        public:
            SCLDecoder(const PolarCodeSpecification* spec, size_t maxPaths);
            ~SCLDecoder() override;

            InfVector Decode(const std::vector<double>& inputLLRs) const override;

        private:
            const PolarCodeSpecification* mSpec;
            size_t mNLayers;
            size_t mMaxPaths;

            mutable std::vector<double> mPathMetrics;

            template <typename T>
            using Vec2d = std::vector<std::vector<T>>;

            mutable Vec2d<double*> mPathLLRs;
            mutable Vec2d<bool*> mPathSymbols;

            mutable Vec2d<size_t> mArrayReferenceCount;
            mutable Vec2d<size_t> mPathToArray;
            mutable Vec2d<size_t> mInactiveArrayIndices;

            mutable std::vector<size_t> mInactivePathIndices;
            mutable std::vector<bool> mIsActivePath;
        
            inline size_t GetLayerSize(size_t layer) const;

            void Initialize() const;

            size_t ClonePath(size_t pathIdx) const;
            void KillPath(size_t pathIdx) const;

            size_t GetArrayIndex(size_t layer, size_t pathIdx) const;

            double* GetLLRsArray(size_t layer, size_t pathIdx) const;
            bool* GetSymbolsArray(size_t layer, size_t pathIdx) const;

            void ProcessLayer(size_t layer, size_t phase) const;
            void ProcessPhase(size_t phase) const;

            void ProcessNegativeBranch(size_t layer, size_t phase) const;
            void ProcessPositiveBranch(size_t layer, size_t phase) const;

            void ContinuePaths_Frozen(size_t phase) const;
            void ContinuePaths_Unfrozen(size_t phase) const;

            void AssignSymbol(size_t pathIdx, bool symbol) const;
        };
    } // namsespace PolarCode
} // namespace Codec
