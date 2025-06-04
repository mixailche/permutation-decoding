#include <iostream>
#include <stdexcept>
#include <fstream>
#include <array>

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

static constexpr size_t Poly(const std::vector<size_t>& degrees)
{
    size_t mask = 0;
    for (size_t deg : degrees) {
        mask |= 1ull << deg;
    }
    return mask;
}

static codec::PolarSpecification BuildPolarSubcode(
    size_t numLayers, size_t dimension,
    const std::vector<double>& errorProbs,
    const running::ArgsReader& reader)
{
    static constexpr std::array<size_t, 15> primitivePolynomials = {
        /* GF(2^2): x^2+x+1 */            Poly({ 2, 1, 0 }),
        /* GF(2^3): x^3+x+1 */            Poly({ 3, 1, 0 }),
        /* GF(2^4): x^4+x+1 */            Poly({ 4, 1, 0 }),
        /* GF(2^5): x^5+x^2+1 */          Poly({ 5, 2, 0 }),
        /* GF(2^6): x^6+x+1 */            Poly({ 6, 1, 0 }),
        /* GF(2^7): x^7+x^3+1 */          Poly({ 7, 3, 0 }),
        /* GF(2^8): x^8+x^4+x^3+x^2+1 */  Poly({ 8, 4, 3, 2, 0 }),
        /* GF(2^9): x^9+x^4+1 */          Poly({ 9, 4, 1 }),
        /* GF(2^10): x^10+x^3+1 */        Poly({ 10, 3, 0 }),
        /* GF(2^11): x^11+x^2+1 */        Poly({ 11, 2, 0 }),
        /* GF(2^12): x^12+x^6+x^4+x+1 */  Poly({ 12, 6, 4, 1, 0 }),
        /* GF(2^13): x^13+x^4+x^3+x+1 */  Poly({ 13, 4, 3, 1, 0 }),
        /* GF(2^14): x^14+x^10+x^6+x+1 */ Poly({ 14, 10, 6, 1, 0 }),
        /* GF(2^15): x^15+x+1 */          Poly({ 15, 1, 0 }),
        /* GF(2^16): x^16+x^12+x^3+x+1 */ Poly({ 16, 12, 3, 1, 0 })
    };

    auto code = reader.GetStringArg("code");

    if (code == "ebch") {
        auto field = math::GField(primitivePolynomials[numLayers - 2]);
        auto minDist = reader.GetNumberArg("d");
        return construct::BuildEBCHSubcode(numLayers, dimension, minDist, field, errorProbs);
    }
    else if (code == "ebch-perm") {
        auto field = math::GField(primitivePolynomials[numLayers - 2]);
        auto minDist = reader.GetNumberArg("d");
        auto blockSizes = reader.GetNumberListArg("blocks");
        return construct::BuildPermFriendlyEBCHSubcode(
            numLayers, dimension, minDist, field, errorProbs, blockSizes
        );
    }
    else if (code == "rand") {
        auto numDFS_A = reader.GetNumberArg("a");
        auto numDFS_B = reader.GetNumberArg("b");
        auto base = construct::BuildFrozenSet(errorProbs, dimension + numDFS_A);
        return construct::BuildRandomizedPolarSubcode(base, errorProbs, numDFS_A, numDFS_B);
    }
    else if (code == "rand-perm") {
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
    else {
        throw std::invalid_argument("Unknown code cunstruction: `" + code + "`");
    }
}

// Usage: Builder.exe <options>
// -out -code -a -b -d -len -dim -BEC -gauss -blocks
// -code: rand, rand-perm, ebch, ebch-perm
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
