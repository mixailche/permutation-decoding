#include <iostream>
#include <fstream>
#include <optional>
#include <string>
#include <unordered_map>
#include <iomanip>

#include "construct/PermSet.h"
#include "construct/FrozenSet.h"
#include "construct/PolarSubcode.h"
#include "codec/SCLDecoder.h"
#include "codec/SCDecoder.h"
#include "codec/PermSCDecoder.h"
#include "codec/PermSCLDecoder.h"
#include "codec/CRCPermSCLDecoder.h"
#include "math/MathUtils.h"
#include "math/MatGF2.h"
#include "math/GField.h"
#include "utils/Utils.hpp"
#include "running/ArgsReader.hpp"
#include "running/Simulate.hpp"

static codec::PolarSpecification ReadSpecification(std::istream& istr)
{
    size_t length, numSymbols;
    istr >> length;
    codec::PolarSpecification spec;

    try {
        spec = codec::PolarSpecification(length);
        istr >> spec.Dimension;
    }
    catch (...) {
        std::cout << "Cannot allocate memory for a code of length " << length << "\n";
        throw;
    }

    while (istr >> numSymbols) {
        std::vector<size_t> equation(numSymbols);
        for (size_t i = 0; i < numSymbols; i++) {
            istr >> equation[i];
        }

        size_t frozenSymbol = equation.back();
        if (numSymbols == 1) {
            spec.StaticFrozen[frozenSymbol] = true;
            continue;
        }

        spec.Dynamic.Frozen[frozenSymbol] = true;
        for (size_t i = 0; i < numSymbols - 1; i++) {
            spec.Dynamic.ForwardEquations[equation[i]].push_back(frozenSymbol);
        }
    }

    return spec;
}

static running::SimulationResult RunSimulator(
    const running::ArgsReader& reader,
    const codec::PolarSpecification& spec)
{
    static auto PrintErrorRate = [](size_t numIterations, size_t numErrors) {
        double fer = static_cast<double>(numErrors) / numIterations;
        std::cout << "steps=" << numIterations
            << ", errors=" << numErrors
            << ", FER=" << std::setprecision(5) << std::scientific << fer
            << "\n";
    };

    auto BuildPerms = [&]() -> std::vector<math::Perm> {
        auto numLayers = utils::IntLog2(spec.Length);
        auto blockSizes = reader.GetNumberListArg("blocks");
        auto count = reader.GetNumberArg("nperms");
        auto minDistUTElems = reader.GetNumberArg("du");
        auto minDistPerm = reader.GetNumberArg("dp");

        return construct::BuildBLTAPermSet(count, numLayers, minDistUTElems, minDistPerm, blockSizes);
    };

    auto signalNoiseRatio = reader.GetDoubleArg("snr");
    auto maxIterations = reader.GetDoubleArg("iter");
    auto maxErrors = reader.GetNumberArg("errors");
    auto callbackPeriod = reader.GetNumberArg("period");
    
    running::SimulationOptions options{
        signalNoiseRatio,
        maxIterations,
        maxErrors,
        callbackPeriod
    };

    auto RunWithDecoder = [&](const codec::Decoder& decoder, const std::vector<bool>& crcGenerator = {}) {
        return running::Simulate(&spec, decoder, options, PrintErrorRate, crcGenerator);
    };

    auto algorithm = reader.GetStringArg("alg");

    if (algorithm == "sc") {
        return RunWithDecoder(codec::SCDecoder(&spec));
    }
    else if (algorithm == "scl") {
        return RunWithDecoder(codec::SCLDecoder(&spec, reader.GetNumberArg("L")));
    }
    else if (algorithm == "perm-sc") {
        return RunWithDecoder(codec::PermSCDecoder(spec, BuildPerms()));
    }
    else if (algorithm == "perm-scl") {
        return RunWithDecoder(codec::PermSCLDecoder(spec, BuildPerms(), reader.GetNumberArg("L")));
    }
    else if (algorithm == "crc-perm-scl") {
        auto rawCrcGenerator = reader.GetNumberListArg("crc");
        std::vector<bool> crcGenerator(rawCrcGenerator.size());
        std::transform(rawCrcGenerator.begin(), rawCrcGenerator.end(), crcGenerator.begin(),
            [](size_t x) -> bool { return x; });
        return RunWithDecoder(codec::CRCPermSCLDecoder(spec, crcGenerator, BuildPerms(), reader.GetNumberArg("L")), crcGenerator);
    }
    else {
        throw std::invalid_argument("Unkown algorithm");
    }
}

// Usage: Decoder.exe <parameters>
// -spec -alg -snr -iter -errors -period -L -du -dp -nperms -blocks
int main(int argc, char** argv)
{
    running::ArgsReader reader(argc, argv, {
        { "iter", "10000000" },
        { "errors", "100" },
        { "period", "100" },
        { "du", "0" },
        { "dp", "0" }
    });

    try {
        codec::PolarSpecification spec;
        if (auto specStrm = std::ifstream(reader.GetStringArg("spec"))) {
            spec = ReadSpecification(specStrm);
        }
        else {
            std::cout << "Cannot open specification file\n";
            return 1;
        }

        auto result = RunSimulator(reader, spec);
        std::cout
            << "FER: " << std::setprecision(5) << std::scientific << result.FrameErrorRate << "\n"
            << "Average number of operations: " << result.AvgNumOperations << "\n";
    }
    catch (const std::invalid_argument& ex) {
        std::cout << ex.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cout << "Unexpected error\n";
        return 1;
    }

    return 0;
}
