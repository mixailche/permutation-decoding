#include <vector>
#include "codec/Decoder.h"
#include "codec/PolarSpecification.h"

namespace codec {
    class SCLDecoder : public Decoder {
    public:
        SCLDecoder(const PolarSpecification* spec, size_t maxPaths);
        ~SCLDecoder() override;

        std::vector<bool> Decode(const std::vector<double>& inputLLRs) const override;
        size_t NumOperations() const override;

    private:
        const PolarSpecification* mSpec;
        size_t mNumLayers;
        size_t mMaxPaths;

        mutable size_t mNumOperations;
        mutable std::vector<double> mPathMetrics;

        template <typename T>
        using Vec2d = std::vector<std::vector<T>>;

        mutable Vec2d<double*> mPathLLRs;
        mutable Vec2d<bool*> mPathSymbols;
        mutable std::vector<bool*> mPathDynFrozenSymbols;

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

        double GetContinuedPathMetric(double currentMetric, double y, bool u) const;

        void ProcessNegativeBranch(size_t layer, size_t phase) const;
        void ProcessPositiveBranch(size_t layer, size_t phase) const;

        void ContinuePaths_StaticFrozen(size_t phase) const;
        void ContinuePaths_DynamicFrozen(size_t phase) const;
        void ContinuePaths_Unfrozen(size_t phase) const;

        void AssignSymbol(size_t pathIdx, size_t phase, bool symbol) const;
    };
}
