#include <numeric>
#include <stdexcept>
#include <optional>
#include <queue>
#include <utility>
#include <random>
#include <set>
#include <unordered_set>
#include <algorithm>

#include "utils/Utils.hpp"
#include "construct/PolarSubcode.h"
#include "math/Mathutils.h"

using math::GField;
using math::MatGF2;
using math::VecGF2;
using math::Perm;

static MatGF2 BuildExtendedFieldCheckMatrix(
    const GField& field,
    const std::vector<size_t>& rootDegrees)
{
    size_t numDigits = field.NumDigits();
    auto checkMatrix = MatGF2::Zeros(numDigits * rootDegrees.size(), field.Size());

    for (size_t rowBlock = 0; rowBlock < rootDegrees.size(); rowBlock++) {
        auto deg = rootDegrees[rowBlock];

        for (size_t col = 0; col < checkMatrix.NumCols(); col++) {
            auto power = field.Pow(col, deg);
            for (size_t i = 0; i < numDigits; i++) {
                checkMatrix.Set(numDigits * rowBlock + i, col, power & (1ull << i));
            }
        }
    }

    return checkMatrix;
}

static MatGF2 BuildEBCHCheckMatrix(const GField& field, size_t minDist)
{
    return BuildExtendedFieldCheckMatrix(field, utils::Iota(minDist - 1));
}

codec::PolarSpecification construct::BuildEBCHSubcode(
    size_t numLayers, size_t dimension, size_t minDist,
    const math::GField& field,
    const std::vector<double>& errorProbs)
{
    auto length = 1ull << numLayers;
    auto checkMatrix = BuildEBCHCheckMatrix(field, minDist);
    auto kernel = math::BuildArikanKernel(numLayers);
    auto freezingMatrix = checkMatrix * kernel.Transpose();
    freezingMatrix.RunGaussianElimination();
    codec::PolarSpecification spec(freezingMatrix);
    
    if (spec.Dimension < dimension) {
        throw std::invalid_argument("There is no EBCH subcode of given length and dimension");
    }

    std::vector<std::vector<size_t>> equations(length);
    for (size_t i = 0; i < length; i++) {
        for (size_t j : spec.Dynamic.ForwardEquations[i]) {
            equations[j].push_back(i);
        }
    }

    auto indices = utils::SortingPerm(errorProbs, std::greater<double>());
    if (spec.Dimension < dimension) {
        throw std::invalid_argument("There is no EBCH subcode of the given parameters");
    }
    auto needToFreeze = spec.Dimension - dimension;

    for (auto i : indices) {
        if (spec.StaticFrozen[i] || spec.Dynamic.Frozen[i]) {
            continue;
        }
        if (needToFreeze-- == 0) {
            break;
        }
        spec.Dimension--;
        spec.Dynamic.ForwardEquations[i].clear();
        spec.StaticFrozen[i] = true;
    }

    for (size_t i = 0; i < length; i++) {
        auto& eq = equations[i];
        if (spec.Dynamic.Frozen[i] && std::all_of(eq.begin(), eq.end(),
            [&](size_t j) { return spec.StaticFrozen[j]; }))
        {
            spec.Dynamic.ForwardEquations[i].clear();
            spec.Dynamic.Frozen[i] = false;
            spec.StaticFrozen[i] = true;
        }
    }

    return spec;
}

codec::PolarSpecification construct::BuildPermFriendlyEBCHSubcode(
    size_t numLayers, size_t dimension, size_t minDist,
    const GField& field,
    const std::vector<double>& errorProbs,
    const std::vector<size_t>& blockSizes)
{
    auto length = 1ull << numLayers;
    auto checkMatrix = BuildEBCHCheckMatrix(field, minDist);
    auto kernel = math::BuildArikanKernel(numLayers);
    auto freezingMatrix = checkMatrix * kernel.Transpose();
    freezingMatrix.RunGaussianElimination();
    codec::PolarSpecification spec(freezingMatrix);

    if (spec.Dimension < dimension) {
        throw std::invalid_argument("There is no EBCH subcode of given length and dimension");
    }

    std::vector<std::vector<size_t>> equations(length);
    for (size_t i = 0; i < length; i++) {
        for (size_t j : spec.Dynamic.ForwardEquations[i]) {
            equations[j].push_back(i);
        }
    }

    auto rawPerms = utils::ListBlockDigitsPermutations(blockSizes);
    std::vector<Perm> perms(rawPerms.size());
    std::transform(rawPerms.begin(), rawPerms.end(), perms.begin(), Perm::MakeDigits);

    std::vector classified(length, false);
    std::vector<std::vector<size_t>> classes;
    
    for (size_t i = 0; i < length; i++) {
        if (classified[i]) {
            continue;
        }
        std::unordered_set<size_t> clazz;
        for (const auto& perm : perms) {
            auto j = perm[i];
            classified[j] = true;
            if (!spec.StaticFrozen[j] && !spec.Dynamic.Frozen[j]) {
                clazz.insert(j);
            }
        }
        classes.emplace_back(clazz.begin(), clazz.end());
    }

    std::vector<size_t> weights(classes.size());
    std::transform(classes.begin(), classes.end(), weights.begin(),
        [](const std::vector<size_t>& clazz) { return clazz.size(); });

    std::vector<double> costs(classes.size());
    std::transform(classes.begin(), classes.end(), costs.begin(),
        [&](const std::vector<size_t>& clazz) {
            double cost = 0;
            for (size_t i : clazz) {
                cost += std::log(1 - errorProbs[i]);
            }
            return cost;
        });

    auto needToFreeze = spec.Dimension - dimension;
    auto solution = utils::SolveKnapsackProblem(classes.size(), needToFreeze, weights, costs);
    for (auto pickedClassIdx : solution) {
        for (size_t i : classes[pickedClassIdx]) {
            spec.Dimension--;
            spec.Dynamic.ForwardEquations[i].clear();
            spec.StaticFrozen[i] = true;
        }
    }

    for (size_t i = 0; i < length; i++) {
        auto& eq = equations[i];
        if (spec.Dynamic.Frozen[i] && std::all_of(eq.begin(), eq.end(),
            [&](size_t j) { return spec.StaticFrozen[j]; }))
        {
            spec.Dynamic.ForwardEquations[i].clear();
            spec.Dynamic.Frozen[i] = false;
            spec.StaticFrozen[i] = true;
        }
    }

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

codec::PolarSpecification construct::BuildRandomizedPolarSubcode(
    const std::vector<bool>& baseCodeFrozenSymbols,
    const std::vector<double>& errorProbs, size_t numDFS_A, size_t numDFS_B)
{
    auto length = baseCodeFrozenSymbols.size();
    auto dimension = std::count(baseCodeFrozenSymbols.begin(), baseCodeFrozenSymbols.end(), false);

    codec::PolarSpecification spec(length);
    spec.Dimension = dimension - numDFS_A;

    std::vector<size_t> unfrozen, frozen;
    for (size_t i = 0; i < length; i++) {
        if (baseCodeFrozenSymbols[i]) {
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
            if (!spec.StaticFrozen[j] && !spec.Dynamic.Frozen[j] && utils::RandomGF2()) {
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
        auto wLeft = utils::HammingWeight(i);
        auto wRight = utils::HammingWeight(j);
        return wLeft < wRight || wLeft == wRight && i > j;
    });
    std::sort(frozen.begin(), frozen.end(), [&](size_t i, size_t j) {
        return errorProbs[i] < errorProbs[j];
    });

    std::for_each_n(unfrozen.begin(), numDFS_A, applyDFC);
    std::for_each_n(frozen.begin(), numDFS_B, applyDFC);

    return spec;
}

codec::PolarSpecification construct::BuildPermFriendlyRandomizedPolarSubcode(
    const std::vector<bool>& baseCodeFrozenSymbols,
    const std::vector<double>& errorProbs,
    const std::vector<utils::RawPerm>& linearPerms,
    size_t numDFS_A, size_t numDFS_B)
{
    auto length = baseCodeFrozenSymbols.size();
    auto dimension = std::count(baseCodeFrozenSymbols.begin(), baseCodeFrozenSymbols.end(), false);

    codec::PolarSpecification spec(length);
    spec.Dimension = dimension - numDFS_A;

    std::vector<size_t> unfrozen, frozen;
    for (size_t i = 0; i < spec.Length; i++) {
        if (baseCodeFrozenSymbols[i]) {
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
            if (!spec.StaticFrozen[j] && !spec.Dynamic.Frozen[j] && utils::RandomGF2()
                && std::all_of(linearPerms.begin(), linearPerms.end(), [&](const auto& perm) {
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
        auto wLeft = utils::HammingWeight(i);
        auto wRight = utils::HammingWeight(j);
        return wLeft < wRight || wLeft == wRight && i > j;
    });
    std::sort(frozen.begin(), frozen.end(), [&](size_t i, size_t j) {
        return errorProbs[i] < errorProbs[j];
    });

    std::for_each_n(unfrozen.begin(), numDFS_A, applyDFC);
    std::for_each_n(frozen.begin(), numDFS_B, applyDFC);

    return spec;
}

codec::PolarSpecification construct::BuildPlotkinPolarSubcode(
    const codec::PolarSpecification& code1,
    const codec::PolarSpecification& code2)
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

    codec::PolarSpecification spec(length);
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
