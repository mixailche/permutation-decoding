#include <vector>
#include <algorithm>
#include <random>

#include "codec/PolarSpecification.h"
#include "codec/PolarEncoder.h"
#include "codec/Decoder.h"
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
        SimulationOptions opt, Callback callback)
    {
        static std::random_device device;
        static std::mt19937 gen(device());
        static std::uniform_int_distribution<int> symbolDistr(0, 1);
        static auto generateSymbol = [] { return symbolDistr(gen); };

        auto stddev = std::pow(10.0, -opt.SignalNoiseRatio / 20);
        auto variance = stddev * stddev;
        std::normal_distribution<double> noiseDistr(0, stddev);

        auto simulateNoise = [&](bool symbol) -> double {
            auto x = (symbol ? 1 : -1) + noiseDistr(gen);
            return -2 * x / variance;
        };

        std::vector<bool> symbols(pSpec->Dimension);
        std::vector<double> received(pSpec->Length);

        codec::PolarEncoder encoder(pSpec);
        SimulationResult result{};
        size_t numErrors = 0, numIterations;

        for (numIterations = 0; numIterations < opt.MaxIterations; numIterations++) {
            std::generate(symbols.begin(), symbols.end(), generateSymbol);
            auto codeword = encoder.Encode(symbols);

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

