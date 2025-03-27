#include <random>
#include <vector>
#include <algorithm>

#include "Common.h"
#include "Specification.h"

#include "Construct/Perm.hpp"

#include "PolarCode/Encoder.h"
#include "PolarCode/SCDecoder.h"
#include "PolarCode/SCLDecoder.h"
#include "PolarCode/PermSCDecoder.h"
#include "PolarCode/BPDecoder.h"
#include "PolarCode/AutBPDecoder.h"
#include "PolarCode/PermBPDecoder.h"
#include "PolarCode/BPPermSCDecoder.h"

#include "PolarSubcode/Encoder.h"
#include "PolarSubcode/SCLDecoder.h"
#include "PolarSubcode/AutSCDecoer.h"
#include "PolarSubcode/PermSCDecoder.h"
#include "PolarSubcode/PermSCLDecoder.h"

namespace Running {
    namespace detail {
        template <typename Callback>
        void Simulate(
            double snr, size_t nIterations, size_t maxErrors,
            size_t length, size_t dimension,
            const Codec::Encoder& encoder,
            const Codec::Decoder& decoder,
            Callback callback
        );
    } // namespace detail

    template <typename Callback>
    void Simulate_PolarCodeSC(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarCodeSpecification& spec,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarCodeSCL(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarCodeSpecification& spec, size_t maxPaths,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarCodePermSC(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarCodeSpecification& spec,
        std::vector<Construct::Perm> permutations,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarCodeBP(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarCodeSpecification& spec, size_t numIterMax,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarCodeAutBP(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarCodeSpecification& spec,
        size_t numIterMax, std::vector<Perm> permutations,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarCodePermBP(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarCodeSpecification& spec,
        size_t numIterMax, std::vector<Perm> permutations,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarCodeBPPermSC(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarCodeSpecification& spec,
        size_t numIterMax, std::vector<Perm> permutations,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarSubcodeSC(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarSubcodeSpecification& spec,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarSubcodeSCL(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarSubcodeSpecification& spec, size_t maxPaths,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarSubcodeAutSC(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarSubcodeSpecification& spec,
        std::vector<Construct::Perm> permuations,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarSubcodePermSC(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarSubcodeSpecification& spec,
        std::vector<Construct::Perm> permuations,
        Callback callback
    );

    template <typename Callback>
    void Simulate_PolarSubcodePermSCL(
        double snr, size_t nIterations, size_t maxErrors,
        const Codec::PolarSubcodeSpecification& spec,
        std::vector<Construct::Perm> permutations, size_t maxPaths,
        Callback callback
    );
} // namespace running

template <typename Callback>
void Running::detail::Simulate(
    double snr, size_t nIterations, size_t maxErrors,
    size_t length, size_t dimension,
    const Codec::Encoder& encoder,
    const Codec::Decoder& decoder,
    Callback callback)
{
    static std::random_device device;
    static std::mt19937 gen(device());
    static std::uniform_int_distribution<int> symbolDistr(0, 1);
    static auto generateSymbol = [] { return symbolDistr(gen); };

    auto stddev = std::pow(10.0, -snr / 20);
    auto variance = stddev * stddev;
    std::normal_distribution<double> noiseDistr(0, stddev);

    auto simulateNoise = [&](bool symbol) -> double {
        auto x = (symbol ? 1 : -1) + noiseDistr(gen);
        return -2 * x / variance;
    };

    std::vector<bool> symbols(dimension);
    std::vector<double> received(length);

    size_t nErrors = 0;

    for (size_t nSteps = 1, nErrors = 0; nSteps <= nIterations; nSteps++) {
        std::generate(symbols.begin(), symbols.end(), generateSymbol);
        auto codeword = encoder.Encode(symbols);

        std::transform(codeword.begin(), codeword.end(), received.begin(), simulateNoise);
        auto decoded = decoder.Decode(received);

        if (decoded != codeword && ++nErrors == maxErrors) {
            std::cout << "FER=" << static_cast<double>(nErrors) / nSteps << std::endl;
            return;
        }

        callback(nSteps, nErrors);
    }

    std::cout << "FER=" << static_cast<double>(nErrors) / nIterations << std::endl;
}

template<typename Callback>
void Running::Simulate_PolarCodeSC(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarCodeSpecification& spec,
    Callback callback
) {
    Codec::PolarCode::Encoder encoder(&spec);
    Codec::PolarCode::SCDecoder decoder(&spec);

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarCodeSCL(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarCodeSpecification& spec, size_t maxPaths,
    Callback callback
) {
    Codec::PolarCode::Encoder encoder(&spec);
    Codec::PolarCode::SCLDecoder decoder(&spec, maxPaths);

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarCodePermSC(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarCodeSpecification& spec,
    std::vector<Construct::Perm> permutations,
    Callback callback
) {
    Codec::PolarCode::Encoder encoder(&spec);
    Codec::PolarCode::PermSCDecoder decoder(&spec, std::move(permutations));

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarCodeBP(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarCodeSpecification& spec, size_t numIterMax,
    Callback callback)
{
    Codec::PolarCode::Encoder encoder(&spec);
    Codec::PolarCode::BPDecoder decoder(&spec, numIterMax);

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarCodeAutBP(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarCodeSpecification& spec,
    size_t numIterMax, std::vector<Perm> permutations,
    Callback callback)
{
    Codec::PolarCode::Encoder encoder(&spec);
    Codec::PolarCode::AutBPDecoder decoder(&spec, numIterMax, std::move(permutations));

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarCodePermBP(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarCodeSpecification& spec, size_t numIterMax,
    std::vector<Perm> permutations,
    Callback callback)
{
    Codec::PolarCode::Encoder encoder(&spec);
    Codec::PolarCode::PermBPDecoder decoder(&spec, numIterMax, std::move(permutations));

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarCodeBPPermSC(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarCodeSpecification& spec, size_t numIterMax,
    std::vector<Perm> permutations,
    Callback callback)
{
    Codec::PolarCode::Encoder encoder(&spec);
    Codec::PolarCode::BPPermSCDecoder decoder(&spec, numIterMax, std::move(permutations));

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarSubcodeSC(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarSubcodeSpecification& spec,
    Callback callback)
{
    Codec::PolarSubcode::Encoder encoder(&spec);
    Codec::PolarSubcode::SCDecoder decoder(&spec);

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarSubcodeSCL(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarSubcodeSpecification& spec, size_t maxPaths,
    Callback callback)
{
    Codec::PolarSubcode::Encoder encoder(&spec);
    Codec::PolarSubcode::SCLDecoder decoder(&spec, maxPaths);
    
    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarSubcodeAutSC(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarSubcodeSpecification& spec,
    std::vector<Construct::Perm> permuations,
    Callback callback)
{
    Codec::PolarSubcode::Encoder encoder(&spec);
    Codec::PolarSubcode::AutSCDecoder decoder(&spec, std::move(permuations));

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarSubcodePermSC(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarSubcodeSpecification& spec,
    std::vector<Construct::Perm> permuations,
    Callback callback)
{
    Codec::PolarSubcode::Encoder encoder(&spec);
    Codec::PolarSubcode::PermSCDecoder decoder(&spec, std::move(permuations));

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}

template<typename Callback>
void Running::Simulate_PolarSubcodePermSCL(
    double snr, size_t nIterations, size_t maxErrors,
    const Codec::PolarSubcodeSpecification& spec,
    std::vector<Construct::Perm> permutations, size_t maxPaths,
    Callback callback)
{
    Codec::PolarSubcode::Encoder encoder(&spec);
    Codec::PolarSubcode::PermSCLDecoder decoder(&spec, std::move(permutations), maxPaths);

    detail::Simulate(
        snr, nIterations, maxErrors,
        spec.Length, spec.Dimension,
        encoder, decoder, callback
    );
}
