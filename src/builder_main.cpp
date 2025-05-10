#include <iostream>
#include <stdexcept>
#include <fstream>

#include "construct/FrozenSet.h"
#include "construct/PolarSubcode.h"
#include "running/ArgsReader.hpp"

static std::ostream& operator<<(std::ostream& ostr, const codec::PolarSpecification& spec)
{
    ostr << spec.Length << " " << spec.Dimension << "\n";

    std::vector<std::vector<size_t>> equations(spec.Length);
    for (size_t i = 0; i < spec.Length; i++) {
        for (size_t j : spec.Dynamic.ForwardEquations[i]) {
            equations[j].push_back(i);
        }
    }

    for (size_t i = 0; i < spec.Length; i++) {
        if (!spec.StaticFrozen[i] && !spec.Dynamic.Frozen[i]) {
            continue;
        }
        ostr << equations[i].size() + 1;
        for (size_t j : equations[i]) {
            ostr << " " << j;
        }
        ostr << " " << i << "\n";
    }

    return ostr;
}

static codec::PolarSpecification BuildPolarSubcode(
    size_t numLayers, size_t dimension,
    const std::vector<double>& errorProbs,
    const running::ArgsReader& reader)
{
    auto code = reader.GetStringArg("code");

    if (code == "ebch") {
        auto field = math::GField(reader.GetNumberArg("mod"));
        auto minDist = reader.GetNumberArg("d");
        return construct::BuildEBCHSubcode(numLayers, dimension, minDist, field, errorProbs);
    }
    if (code == "rand") {
        auto numDFS_A = reader.GetNumberArg("a");
        auto numDFS_B = reader.GetNumberArg("b");
        auto base = construct::BuildFrozenSet(errorProbs, dimension + numDFS_A);
        return construct::BuildRandomizedPolarSubcode(base, errorProbs, numDFS_A, numDFS_B);
    }
    if (code == "rand-perm") {
        auto numDFS_A = reader.GetNumberArg("a");
        auto numDFS_B = reader.GetNumberArg("b");
        auto blockSizes = reader.GetNumberListArg("blocks");
        auto base = construct::BuildFrozenSet_AutomorphismInvariant(
            numLayers, dimension + numDFS_A, blockSizes, errorProbs
        );
        auto linearPerms = utils::ListBlockDigitsPermutations(blockSizes);
        return construct::BuildPermFriendlyRandomizedPolarSubcode(
            base, errorProbs, linearPerms, numDFS_A, numDFS_B
        );
    }

    throw std::invalid_argument("Unknown code cunstruction: `" + code + "`");
}

// Usage: Builder.exe <options>
// -out -code -a -b -d -len -dim -BEC -gauss -blocks -mod
// -code: rand, rand-perm, ebch
int main(int argc, char** argv)
{
    running::ArgsReader reader(argc, argv);

    try {
        auto length = reader.GetNumberArg("len");
        auto dimension = reader.GetNumberArg("dim");
        auto numLayers = utils::IntLog2(length);
        if (length < dimension) {
            std::cout << "Error: len < dim\n";
        }
        if (length != (1ull << numLayers)) {
            std::cout << "Code length must be a power of 2\n";
            return 1;
        }

        auto isBEC = reader.HasArg("BEC");
        auto isGauss = reader.HasArg("gauss");
        if (isBEC && isGauss) {
            std::cout << "Both -BEC and -gauss arguments specified\n";
            return 1;
        }
        if (!isBEC && !isGauss) {
            std::cout << "Target channel not specified\n";
            return 1;
        }

        auto errorProbs = isBEC
            ? construct::BhattacharyyaParameters(numLayers, reader.GetDoubleArg("BEC"))
            : construct::DensityEvolutionGA(numLayers, reader.GetDoubleArg("gauss"));

        auto spec = BuildPolarSubcode(numLayers, dimension, errorProbs, reader);
        
        std::ofstream fout(reader.GetStringArg("out"));
        fout << spec;
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
