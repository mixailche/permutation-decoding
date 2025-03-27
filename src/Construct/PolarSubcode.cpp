#include <numeric>
#include <stdexcept>
#include <optional>
#include <queue>
#include <utility>
#include <random>
#include <set>
#include <unordered_set>
#include <algorithm>

#include "NTL/GF2E.h"
#include "NTL/mat_GF2.h"
#include "NTLHelp.h"
#include "Utils.hpp"
#include "Construct/PolarSubcode.h"
#include "Construct/Perm.hpp"
#include "Construct/Frozen.h"

using NTL::GF2E;
using NTL::GF2X;
using NTL::GF2;

static NTL::Mat<GF2> BuildExtendedFieldCheckMatrix(const std::vector<size_t> rootDegrees)
{
    size_t numDigits = GF2E::degree();

    NTL::mat_GF2 checkMatrix(NTL::INIT_SIZE,
        numDigits * rootDegrees.size(),
        NTL::conv<long>(GF2E::cardinality()));

    for (size_t rowBlock = 0; rowBlock < rootDegrees.size(); rowBlock++) {
        auto deg = NTL::ZZ(rootDegrees[rowBlock]);
        
        for (size_t col = 0; col < checkMatrix.NumCols(); col++) {
            auto elem = Math::IndexToGF2E(col, numDigits);
            auto power = NTL::conv<GF2X>(NTL::power(elem, deg));

            for (size_t i = 0; i < numDigits; i++) {
                checkMatrix.put(numDigits * rowBlock + i, col, NTL::coeff(power, i));
            }
        }
    }

    return checkMatrix;
}

static NTL::Mat<GF2> BuildEbchCheckMatrix(size_t minDist)
{
    return BuildExtendedFieldCheckMatrix(Utils::Iota(minDist - 1));
}

NTL::mat_GF2 Construct::BuildFreezingMatrix(const Codec::PolarSubcodeSpecification& spec)
{
    NTL::mat_GF2 freezingMatrix(NTL::INIT_SIZE, spec.Length - spec.Dimension, spec.Length);
    std::vector<size_t> dynamicFrozenSymbolRows(spec.Length);

    for (size_t i = 0, j = 0; j < spec.Length; i++, j++) {
        if (spec.StaticFrozen[j]) {
            freezingMatrix.put(i, j, 1);
        }
        else if (spec.Dynamic.Frozen[j]) {
            freezingMatrix.put(i, j, 1);
            dynamicFrozenSymbolRows[j] = i;
        }
        else {
            i--;
        }
    }

    for (size_t i = 0; i < spec.Length; i++) {
        for (auto dyn : spec.Dynamic.ForwardEquations[i]) {
            freezingMatrix[dynamicFrozenSymbolRows[dyn]][i] = 1;
        }
    }

    return freezingMatrix;
}

std::vector<std::vector<size_t>> Construct::BuildDynamicFreezingEquations(
    const NTL::mat_GF2& freezingMatrix,
    std::vector<bool>& staticFrozen)
{
    std::vector<std::vector<size_t>> equations(staticFrozen.size());

    for (size_t r = 0; r < freezingMatrix.NumRows(); r++) {
        auto& row = freezingMatrix[r];
        if (auto end = Math::HighestOnePos(row)) {
            auto i = *end;
            auto& equation = equations[i];

            for (size_t j = 0; j < i; j++) {
                if (NTL::IsOne(row[j])) {
                    equation.push_back(j);
                }
            }

            if (equation.empty()) {
                staticFrozen[i] = true;
            }
        }
    }

    return equations;
}

Codec::PolarSubcodeSpecification Construct::BuildEBCHSubcode(size_t numLayers, size_t dimension, size_t minDist, const std::vector<double>& errorProbs)
{
    auto length = 1ull << numLayers;
    auto needToFreeze = length - dimension;

    auto checkMatrix = BuildEbchCheckMatrix(minDist);
    auto kernel = NTL::transpose(Math::BuildArikanKernel(numLayers));
    auto freezingMatrix = checkMatrix * kernel;
    auto numDynamicFrozen = Math::RunDoubleGaussianElimination(freezingMatrix);

    if (needToFreeze < numDynamicFrozen) {
        throw std::invalid_argument("There is not EBCH subcode with the given parameters");
    }
    needToFreeze -= numDynamicFrozen;

    Codec::PolarSubcodeSpecification spec(length);
    spec.Dimension = dimension;

    auto equations = BuildDynamicFreezingEquations(freezingMatrix, spec.StaticFrozen);
    auto indices = Utils::SortingPerm(errorProbs, std::greater<double>());
    
    for (auto i : indices) {
        if (spec.StaticFrozen[i] || !equations[i].empty()) {
            continue;
        }
        if (needToFreeze-- == 0) {
            break;
        }
        spec.StaticFrozen[i] = true;
    }

    for (size_t i = 0; i < length; i++) {
        for (size_t j : equations[i]) {
            if (!spec.StaticFrozen[j]) {
                spec.Dynamic.Frozen[i] = true;
                spec.Dynamic.ForwardEquations[j].push_back(i);
            }
        }

        if (!equations[i].empty() && !spec.Dynamic.Frozen[i]) {
            // all summands are trivial
            spec.StaticFrozen[i] = true;
        }
    }

    return spec;
}

Codec::PolarSubcodeSpecification Construct::BuildAutFriendlyEBCHSubcode(
    size_t numLayers, size_t dimension, size_t minDist,
    const std::vector<double>& errorProbs,
    const std::vector<size_t>& blockSizes)
{
    static std::random_device device;
    static std::mt19937 gen(device());

    auto length = 1ull << numLayers;
    auto needToFreeze = length - dimension;

    auto checkMatrix = BuildEbchCheckMatrix(minDist);
    auto kernel = NTL::transpose(Math::BuildArikanKernel(numLayers));
    auto freezingMatrix = checkMatrix * kernel;
    auto numDynamicFrozen = Math::RunDoubleGaussianElimination(freezingMatrix);

    if (needToFreeze < numDynamicFrozen) {
        throw std::invalid_argument("There is not EBCH subcode with the given parameters");
    }
    needToFreeze -= numDynamicFrozen;

    Codec::PolarSubcodeSpecification spec(length);
    spec.Dimension = length - numDynamicFrozen;
    auto equations = Construct::BuildDynamicFreezingEquations(freezingMatrix, spec.StaticFrozen);
    
    auto digitPerms = Utils::ListBlockDigitsPermutations(blockSizes);
    std::vector<std::vector<size_t>> classes;
    std::vector<size_t> weights;
    std::vector<double> costs;
    std::vector classified(length, false);

    for (size_t i = 0; i < length; i++) {
        if (classified[i]) {
            continue;
        }
        std::vector<size_t> clazz;
        double cost = 0;
        for (auto& digitsPerm : digitPerms) {
            auto perm = Perm::MakeDigits(digitsPerm);
            auto j = perm[i];
            classified[j] = true;
            if (spec.StaticFrozen[j] || !equations[j].empty()) {
                continue;
            }
            cost += std::log(1 - errorProbs[j]);
            clazz.push_back(j);
        }
        costs.push_back(cost);
        weights.push_back(clazz.size());
        classes.push_back(std::move(clazz));
    }

    auto solution = Utils::SolveKnapsackProblem(classes.size(), needToFreeze, weights, costs);
    for (auto pickedClassIdx : solution) {
        for (auto i : classes[pickedClassIdx]) {
            spec.StaticFrozen[i] = true;
            spec.Dimension--;
        }
    }

    for (size_t i = 0; i < length; i++) {
        for (size_t j : equations[i]) {
            if (!spec.StaticFrozen[j]) {
                spec.Dynamic.Frozen[i] = true;
                spec.Dynamic.ForwardEquations[j].push_back(i);
            }
        }

        if (!equations[i].empty() && !spec.Dynamic.Frozen[i]) {
            // all summands are trivial
            spec.StaticFrozen[i] = true;
        }
    }

    return Codec::PolarSubcodeSpecification(freezingMatrix);
}

// Deprecated
static Codec::PolarSubcodeSpecification BROKEN_BuildCyclicPolarSubcode(
    size_t numLayers, size_t dimension,
    const std::vector<double>& errorProbs)
{
    auto length = 1ull << numLayers;
    std::vector<NTL::mat_GF2> freezingSubmatrices;
    std::vector<double> costs;
    std::vector<size_t> weights;
    std::vector used(length, false);

    for (size_t i = 0; i < length - 1; i++) {
        if (used[i]) {
            continue;
        }
        used[i] = true;
        
        size_t currentElem = i;
        std::vector<size_t> coset = { currentElem };
        do {
            currentElem = (currentElem * 2) % (length - 1);
            coset.push_back(currentElem);
            used[currentElem] = true;
        } while (currentElem != coset.front());

        auto checkMatrix = BuildExtendedFieldCheckMatrix(coset);
        auto kernel = NTL::transpose(Math::BuildArikanKernel(numLayers));
        auto freezingMatrix = checkMatrix * kernel;
        Math::RunDoubleGaussianElimination(freezingMatrix);
        auto staticFrozen = std::vector(length, false);
        auto equations = Construct::BuildDynamicFreezingEquations(freezingMatrix, staticFrozen);
        
        std::cout << "Coset leader: " << i << std::endl;

        double cost = 0;
        for (size_t j = 0; j < length; j++) {
            if (staticFrozen[j] || !equations[j].empty()) {
                cost += std::log(1 - errorProbs[j]);
                std::cout << j << std::endl;
            }
        }
        costs.push_back(cost);
        
        NTL::vec_vec_GF2 nonZeroRows;
        for (size_t j = 0; j < freezingMatrix.NumRows(); j++) {
            auto& row = freezingMatrix[j];
            if (std::any_of(row.begin(), row.end(), [](GF2 x) { return NTL::IsOne(x); })) {
                nonZeroRows.append(row);
            }
        }
        weights.push_back(nonZeroRows.length());
        freezingSubmatrices.push_back(NTL::to_mat_GF2(nonZeroRows));

        std::cout << freezingSubmatrices.back() << std::endl;
    }

    auto solution = Utils::SolveKnapsackProblem(weights.size(), length - dimension, weights, costs);
    auto freezingMatrixRows = NTL::vec_vec_GF2();
    for (auto pickedCosetIdx : solution) {
        for (size_t i = 0; i < freezingSubmatrices[pickedCosetIdx].NumRows(); i++) {
            freezingMatrixRows.append(freezingSubmatrices[pickedCosetIdx][i]);
        }
    }
    auto freezingMatrix = NTL::to_mat_GF2(freezingMatrixRows);
    Math::RunDoubleGaussianElimination(freezingMatrix);

    Codec::PolarSubcodeSpecification spec(length);
    auto equations = Construct::BuildDynamicFreezingEquations(freezingMatrix, spec.StaticFrozen);
    for (size_t i = 0; i < length; i++) {
        for (auto j : equations[i]) {
            spec.Dynamic.Frozen[i] = true;
            spec.Dynamic.ForwardEquations[j].push_back(i);
        }
    }

    spec.Dimension =
        std::count(spec.StaticFrozen.begin(), spec.StaticFrozen.end(), false) -
        std::count(spec.Dynamic.Frozen.begin(), spec.Dynamic.Frozen.end(), true);

    return spec;
}

Codec::PolarSubcodeSpecification Construct::BuildCyclicPolarSubcode(
    size_t numLayers, size_t dimension,
    const std::vector<double>& errorProbs)
{
    auto length = 1ull << numLayers;
    std::vector<NTL::mat_GF2> freezingSubmatrices;
    std::vector<double> costs;
    std::vector<size_t> weights;
    std::vector usedElem(length, false);

    for (size_t i = 0; i < length - 1; i++) {
        if (usedElem[i]) {
            continue;
        }
        usedElem[i] = true;

        size_t currentElem = i;
        std::vector<size_t> coset = { currentElem };
        do {
            currentElem = (currentElem * 2) % (length - 1);
            coset.push_back(currentElem);
            usedElem[currentElem] = true;
        } while (currentElem != coset.front());

        auto checkMatrix = BuildExtendedFieldCheckMatrix(coset);
        auto kernel = NTL::transpose(Math::BuildArikanKernel(numLayers));
        auto freezingMatrix = checkMatrix * kernel;
        Math::RunDoubleGaussianElimination(freezingMatrix);

        NTL::vec_vec_GF2 nonZeroRows;
        for (size_t j = 0; j < freezingMatrix.NumRows(); j++) {
            auto& row = freezingMatrix[j];
            if (std::any_of(row.begin(), row.end(), [](GF2 x) { return NTL::IsOne(x); })) {
                nonZeroRows.append(row);
            }
        }
        freezingSubmatrices.push_back(NTL::to_mat_GF2(nonZeroRows));
    }

    NTL::mat_GF2 freezingMatrix(NTL::INIT_SIZE, 0, length);
    std::vector usedCoset(freezingSubmatrices.size(), false);

    while (freezingMatrix.NumRows() < length - dimension) {
        auto currentNumRows = freezingMatrix.NumRows();
        double bestCost = std::numeric_limits<double>::max();
        NTL::mat_GF2 newFreezingMatrix;
        size_t pickedCosetIdx = -1;

        for (size_t i = 0; i < freezingSubmatrices.size(); i++) {
            if (usedCoset[i]) {
                continue;
            }

            auto& submatrix = freezingSubmatrices[i];
            NTL::mat_GF2 forkMatrix = freezingMatrix;

            forkMatrix.SetDims(currentNumRows + submatrix.NumRows(), length);
            for (size_t r = 0; r < submatrix.NumRows(); r++) {
                for (size_t c = 0; c < length; c++) {
                    forkMatrix[currentNumRows + r][c] = submatrix[r][c];
                }
            }

            double cost = 0;
            Math::RunDoubleGaussianElimination(forkMatrix);
            for (size_t r = 0; r < forkMatrix.NumRows(); r++) {
                if (auto end = Math::HighestOnePos(forkMatrix[r])) {
                    cost += std::log(1 - errorProbs[*end]);
                }
            }

            if (cost < bestCost) {
                bestCost = cost;
                newFreezingMatrix = forkMatrix;
                pickedCosetIdx = i;
            }
        }

        usedCoset[pickedCosetIdx] = true;
        freezingMatrix = std::move(newFreezingMatrix);
    }

    Codec::PolarSubcodeSpecification spec(length);
    auto equations = BuildDynamicFreezingEquations(freezingMatrix, spec.StaticFrozen);
    for (size_t i = 0; i < length; i++) {
        for (auto j : equations[i]) {
            spec.Dynamic.Frozen[i] = true;
            spec.Dynamic.ForwardEquations[j].push_back(i);
        }
    }
    spec.Dimension = length - freezingMatrix.NumRows();

    return spec;
}

static size_t ApplyDigitsPerm(const std::vector<size_t>& perm, size_t number)
{
    size_t result = 0;
    for (size_t i = 0; i < perm.size(); i++) {
        if (number & (1ull << perm[i])) {
            result |= (1ull << i);
        }
    }
    return result;
}

Codec::PolarSubcodeSpecification Construct::BuildRandomizedPolarSubcode(
    const Codec::PolarCodeSpecification& polarSpec,
    const std::vector<double>& errorProbs, size_t numDFS_A, size_t numDFS_B)
{
    static std::random_device device;
    static std::mt19937 gen(device());
    static std::uniform_int_distribution distr(0, 1);

    if (polarSpec.Dimension <= numDFS_A) {
        throw std::invalid_argument("Cannot apply type-A DFC");
    }

    Codec::PolarSubcodeSpecification spec(polarSpec.Length);
    spec.Dimension = polarSpec.Dimension - numDFS_A;

    std::vector<size_t> unfrozen, frozen;
    for (size_t i = 0; i < spec.Length; i++) {
        if (polarSpec.Frozen[i]) {
            frozen.push_back(i);
            spec.StaticFrozen[i] = true;
        }
        else {
            unfrozen.push_back(i);
        }
    }

    auto applyDFC = [&](size_t i) -> void {
        std::vector<size_t> equation;

        for (size_t j = 0; j < i; j++) {
            if (!polarSpec.Frozen[j] && distr(gen)) {
                equation.push_back(j);
            }
        }

        if (!(spec.StaticFrozen[i] = equation.empty())) {
            spec.Dynamic.Frozen[i] = true;
            for (auto j : equation) {
                spec.Dynamic.ForwardEquations[j].push_back(i);
            }
        }
    };

    std::sort(unfrozen.begin(), unfrozen.end(), [](size_t i, size_t j) {
        auto wLeft = Utils::HammingWeight(i);
        auto wRight = Utils::HammingWeight(j);
        return wLeft < wRight || wLeft == wRight && i > j;
    });
    std::sort(frozen.begin(), frozen.end(), [&](size_t i, size_t j) {
        return errorProbs[i] < errorProbs[j];
    });

    std::for_each_n(unfrozen.begin(), numDFS_A, applyDFC);
    std::for_each_n(frozen.begin(), numDFS_B, applyDFC);

    return spec;
}

Codec::PolarSubcodeSpecification Construct::BuildPermFriendlyRandomizedPolarSubcode(
    const Codec::PolarCodeSpecification& polarSpec,
    const std::vector<double>& errorProbs,
    const std::vector<Utils::RawPerm>& linearPerms,
    size_t numDFS_A, size_t numDFS_B, size_t l)
{
    static std::random_device device;
    static std::mt19937 gen(device());
    static std::uniform_int_distribution distr(0, 1);

    if (polarSpec.Dimension <= numDFS_A) {
        throw std::invalid_argument("Cannot apply type-A DFC");
    }

    Codec::PolarSubcodeSpecification spec(polarSpec.Length);
    spec.Dimension = polarSpec.Dimension - numDFS_A;

    std::vector<size_t> unfrozen, frozen;
    for (size_t i = 0; i < spec.Length; i++) {
        if (polarSpec.Frozen[i]) {
            frozen.push_back(i);
            spec.StaticFrozen[i] = true;
        }
        else {
            unfrozen.push_back(i);
        }
    }

    auto applyDFC = [&, mask = (1ull << l) - 1](size_t i) -> void {
        std::vector<size_t> equation;

        for (size_t j = 0; j < i; j++) {
            if ((j | mask) < i &&
                !polarSpec.Frozen[j] && distr(gen) &&
                std::all_of(linearPerms.begin(), linearPerms.end(), [&](const auto& perm) {
                    return ApplyDigitsPerm(perm, j) < ApplyDigitsPerm(perm, i);
                })) {
                equation.push_back(j);
            }
        }

        if (!(spec.StaticFrozen[i] = equation.empty())) {
            spec.Dynamic.Frozen[i] = true;
            for (auto j : equation) {
                spec.Dynamic.ForwardEquations[j].push_back(i);
            }
        }
    };

    std::sort(unfrozen.begin(), unfrozen.end(), [](size_t i, size_t j) {
        auto wLeft = Utils::HammingWeight(i);
        auto wRight = Utils::HammingWeight(j);
        return wLeft < wRight || wLeft == wRight && i > j;
    });
    std::sort(frozen.begin(), frozen.end(), [&](size_t i, size_t j) {
        return errorProbs[i] < errorProbs[j];
    });

    std::for_each_n(unfrozen.begin(), numDFS_A, applyDFC);
    std::for_each_n(frozen.begin(), numDFS_B, applyDFC);

    return spec;
}

Codec::PolarSubcodeSpecification Construct::BuildPlotkinPolarSubcode(
    const Codec::PolarSubcodeSpecification& code1,
    const Codec::PolarSubcodeSpecification& code2)
{
    if (code1.Length != code2.Length) {
        throw std::invalid_argument("Codes must be of an equal length");
    }

    auto half = code1.Length;
    auto length = half * 2;

    std::vector<std::vector<size_t>> equations1(half);
    std::vector<std::vector<size_t>> equations2(half);
    for (size_t i = 0; i < half; i++) {
        for (auto j : code1.Dynamic.ForwardEquations[i]) {
            equations1[j].push_back(i);
        }
        for (auto j : code2.Dynamic.ForwardEquations[i]) {
            equations2[j].push_back(i);
        }
    }

    Codec::PolarSubcodeSpecification spec(length);
    spec.Dimension = code1.Dimension + code2.Dimension;

    for (size_t i = 0; i < half; i++) {
        for (auto j : equations1[i]) {
            spec.Dynamic.Frozen[2 * i + 1] = true;
            spec.Dynamic.ForwardEquations[2 * j + 1].push_back(2 * i + 1);
        }
        for (auto j : equations2[i]) {
            spec.Dynamic.Frozen[2 * i] = true;
            spec.Dynamic.ForwardEquations[2 * j].push_back(2 * i);
        }

        if (code1.StaticFrozen[i] && code2.StaticFrozen[i]) {
            spec.StaticFrozen[2 * i] = spec.StaticFrozen[2 * i + 1] = true;
        }
        else if (code1.StaticFrozen[i] || code1.Dynamic.Frozen[i]) {
            spec.Dynamic.Frozen[2 * i + 1] = true;
            spec.Dynamic.ForwardEquations[2 * i].push_back(2 * i + 1);
        }
        else if (code2.StaticFrozen[i]) {
            spec.StaticFrozen[2 * i] = true;
        }
    }

    return spec;
}
