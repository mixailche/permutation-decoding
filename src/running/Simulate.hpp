#include <vector>
#include <algorithm>
#include <random>

#include "codec/PolarSpecification.h"
#include "codec/PolarEncoder.h"
#include "codec/Decoder.h"
#include "math/MathUtils.h"
#include "utils/Utils.hpp"

namespace running {
    struct SimulationOptions {
        double SignalNoiseRatio;
        size_t MaxIterations;
        size_t MaxErrors;
        size_t CallbackPeriod;
    };

    struct SimulationResult {
        double FrameErrorRate;
        double AvgNumOperations;
    };

    // Callback is a function like
    // void DoSomething(size_t numIterations, size_t numErrors);
    // It is being called every opt.CallbackPeriod iterations
    template <typename Callback>
    SimulationResult Simulate(
        const codec::PolarSpecification* pSpec, const codec::Decoder& decoder,
        SimulationOptions opt, Callback callback, const std::vector<bool> crcGenerator = {})
    {
        static std::random_device device;
        static std::mt19937 gen(device());
        static std::uniform_int_distribution<int> symbolDistr(0, 1);
        static auto generateSymbol = [] { return symbolDistr(gen); };

        auto numCRCBits = crcGenerator.empty() ? 0 : crcGenerator.size() - 1;
        auto rate = static_cast<double>(pSpec->Dimension - numCRCBits) / pSpec->Length;
        auto stddev = std::pow(10.0, -opt.SignalNoiseRatio / 20) / std::sqrt(2 * rate);
        auto variance = stddev * stddev;
        std::normal_distribution<double> noiseDistr(0, stddev);

        auto simulateNoise = [&](bool symbol) -> double {
            auto x = (symbol ? 1 : -1) + noiseDistr(gen);
            return -2 * x / variance;
        };

        std::vector<bool> infVector;
        std::vector<double> received(pSpec->Length);

        codec::PolarEncoder encoder(pSpec);
        SimulationResult result{};
        size_t numErrors = 0, numIterations;

        for (numIterations = 0; numIterations < opt.MaxIterations; numIterations++) {
            infVector.resize(pSpec->Dimension - numCRCBits);
            std::generate(infVector.begin(), infVector.end(), generateSymbol);

            std::vector<bool> crc;
            if (!crcGenerator.empty()) {
                crc = math::CalculateCRC(infVector, crcGenerator);
            }

            infVector.insert(infVector.end(), crc.begin(), crc.end());
            auto codeword = encoder.Encode(infVector);

            std::transform(codeword.begin(), codeword.end(), received.begin(), simulateNoise);
            auto decoded = decoder.Decode(received);

            result.AvgNumOperations =
                (result.AvgNumOperations * numIterations + decoder.NumOperations()) / (numIterations + 1);

            if (decoded != codeword && ++numErrors == opt.MaxErrors) {
                break;
            }

            if ((numIterations + 1) % opt.CallbackPeriod == 0) {
                callback(numIterations + 1, numErrors);
            }
        }

        result.FrameErrorRate = static_cast<double>(numErrors) / (numIterations + 1);
        return result;
    }
}

