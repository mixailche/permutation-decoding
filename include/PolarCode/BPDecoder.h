#pragma once

#include <utility>
#include <functional>

#include "PolarCode/Encoder.h"
#include "Common.h"
#include "Specification.h"

namespace Codec {
    namespace PolarCode {
        class BPDecoder : public Decoder {
        public:
            BPDecoder(const PolarCodeSpecification* spec, size_t numIterMax);
            BPDecoder(const PolarCodeSpecification* spec, size_t numIterMax, std::vector<size_t> schedule);
            ~BPDecoder() override;

            InfVector Decode(const std::vector<double>& inputLLRs) const override;
            
            std::pair<bool, std::vector<double>>
            SoftDecode(const std::vector<double>& inputLLRs) const;

        private:
            const PolarCodeSpecification* mSpec;
            size_t mNumLayers;
            size_t mNumIterMax;
            std::vector<size_t> mSchedule;
            PolarCode::Encoder mEncoder;

            mutable std::vector<double*> mLeftMessages;
            mutable std::vector<double*> mRightMessages;

            inline size_t GetLayerSize(size_t layer) const;

            void Initialize(const std::vector<double>& inputLLRs) const;

            void GetInputMessages(
                size_t leftLayer, size_t rightLayer, size_t i1, size_t i2,
                double& l1, double& l2, double& r1, double& r2
            ) const;

            std::vector<std::pair<size_t, size_t>> GetSymbolNodes(size_t layer) const;
            void RunBPIteration() const;

            bool CheckStopCondition() const;
            bool CheckStopCondition(/* output */ Codec::InfVector& codeword) const;

            template <typename T>
            inline std::vector<T> TransformLayer(
                size_t layer, std::function<T(double, double)> func
            ) const;

            inline std::vector<double> GetLayerLLRs(size_t layer) const;
            inline Codec::InfVector GetLayerDecisions(size_t layer) const;
            inline Codec::InfVector GetInfVector() const;
            inline Codec::InfVector GetCodeword() const;
        };
    } // namespace PolarCode
} // namespace Codec
