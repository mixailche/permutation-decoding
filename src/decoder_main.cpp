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

#include "math/MathUtils.h"
#include "math/MatGF2.h"
#include "math/GField.h"
#include "utils/Utils.hpp"
#include "running/ArgsReader.hpp"
#include "running/Simulate.hpp"

class ArgsReader1 {
public:
    ArgsReader1(int argc, char** argv);

    std::optional<std::string> Read();
    std::optional<std::pair<std::string, std::string>> ReadNamed();

private:
    const size_t mNumArgs;
    char** mArgs;
    size_t mIndex;
};

ArgsReader1::ArgsReader1(int argc, char** argv)
    : mNumArgs(argc)
    , mArgs(argv)
    , mIndex(1)
{}

std::optional<std::string> ArgsReader1::Read()
{
    if (mIndex == mNumArgs) {
        return std::nullopt;
    }
    return mArgs[mIndex++];
}

std::optional<std::pair<std::string, std::string>> ArgsReader1::ReadNamed()
{
    auto name = Read();
    auto value = Read();

    if (name && value && !name->empty() && name->at(0) == '-') {
        return std::make_pair(name->substr(1), value.value());
    }
    return std::nullopt;
}

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

    size_t numFrozen;
    for (numFrozen = 0; istr >> numSymbols; numFrozen++) {
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

    auto RunWithDecoder = [&](const codec::Decoder& decoder) {
        return running::Simulate(&spec, decoder, options, PrintErrorRate);
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
    else {
        throw std::invalid_argument("Unkown algorithm");
    }
}

static void PrintWeightsDistr(const codec::PolarSpecification& spec)
{
    auto numLayers = utils::IntLog2(spec.Length);
    auto length = spec.Length;

    std::vector<size_t> weights(numLayers + 1);

    for (size_t i = 0; i < length; i++) {
        if (spec.StaticFrozen[i] || spec.Dynamic.Frozen[i]) {
            weights[math::IndexToVecGF2(i, numLayers).HammingWeight()]++;
        }
    }

    for (size_t w = 0; w <= numLayers; w++) {
        std::cout << w << " -> " << weights[w] << "\n";
    }

    std::cout << std::endl;
}

// Usage: Decoder.exe <parameters>
// -spec -alg -snr -iter -errors -periods -L -du -dp -nperms

int main(int argc, char** argv)
{
    running::ArgsReader reader(argc, argv, {
        { "iter", "1000000" },
        { "errors", "100" },
        { "period", "100" },
        { "du", "0" },
        { "dp", "0" }
    });

    try {
        codec::PolarSpecification spec;
        {
            std::ifstream specStrm(reader.GetStringArg("spec"));
            spec = ReadSpecification(specStrm);
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
