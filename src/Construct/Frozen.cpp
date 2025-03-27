#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <unordered_set>

#include "Construct/Frozen.h"
#include "Construct/Perm.hpp"
#include "Utils.hpp"

std::vector<bool> Construct::BuildFrozenSet(const std::vector<double>& errorProbs, size_t dimension)
{
    auto indices = Utils::SortingPerm(errorProbs, std::greater<double>());

    auto length = errorProbs.size();
    auto need = length - dimension;

    std::vector frozen(length, false);
    for (size_t i = 0; need-- > 0; i++) {
        frozen[indices[i]] = true;
    }

    return frozen;
}

std::vector<bool> Construct::BuildFrozenSet_Bhattacharyya(
    double p, size_t nLayers, size_t dimension)
{
    return BuildFrozenSet(BhattacharyyaParameters(p, nLayers), dimension);
}

std::vector<bool> Construct::BuildFrozenSet_GA(
    double snr, size_t nLayers, size_t dimension)
{
    return BuildFrozenSet(DensityEvolutionGA(snr, nLayers), dimension);
}

static std::vector<bool> BuildFrozenSet_KnapsackProblemSolution(
    size_t length, size_t dimension,
    const std::vector<std::vector<size_t>>& classes,
    const std::vector<double>& errorProbs)
{
    auto numClasses = classes.size();

    std::vector<double> costs(numClasses, 0);
    std::vector<size_t> weights(numClasses);

    for (size_t i = 0; i < numClasses; i++) {
        auto& clazz = classes[i];
        weights[i] = clazz.size();
        for (auto j : clazz) {
            costs[i] += std::log(1 - errorProbs[j]);
        }
    }

    auto solution = Utils::SolveKnapsackProblem(numClasses, length - dimension, weights, costs);
    std::vector frozen(length, false);

    for (auto pickedClassIdx : solution) {
        for (auto i : classes[pickedClassIdx]) {
            frozen[i] = true;
        }
    }

    return frozen;
}

std::vector<bool> Construct::BuildFrozenSet_BlockPermInvariant(
    size_t numLayer, size_t dimension, size_t blockSize,
    const std::vector<double>& errorProbs)
{
    auto length = 1ull << numLayer;
    std::vector<std::vector<size_t>> classes;

    for (size_t classBegin = 0; classBegin < length; classBegin += blockSize) {
        std::vector<size_t> clazz;
        for (size_t i = 0; i < blockSize; i++) {
            if (classBegin + i == length) {
                break;
            }
            clazz.push_back(classBegin + i);
        }
        classes.push_back(std::move(clazz));
    }

    return BuildFrozenSet_KnapsackProblemSolution(length, dimension, classes, errorProbs);
}

template<typename T, typename EvenFunc, typename OddFunc>
static std::vector<T> RunLayersDescent(size_t nLayers, T start, EvenFunc even, OddFunc odd)
{
    std::vector currentLayer = { start };

    for (size_t l = 0; l < nLayers; l++)
    {
        auto half = 1ull << l;
        std::vector<T> nextLayer(half * 2);
        for (size_t i = 0; i < half; i++)
        {
            auto value = currentLayer[i];
            nextLayer[2 * i] = even(value);
            nextLayer[2 * i + 1] = odd(value);
        }
        currentLayer = std::move(nextLayer);
    }

    return currentLayer;
}

std::vector<double> Construct::BhattacharyyaParameters(double p, size_t nLayers)
{
    return RunLayersDescent(nLayers, p,
        [](double x) { return x * (2 - x); },
        [](double x) { return x * x; });
}

static double GaussianCdf(double x)
{
    // https://www.johndcook.com/cpp_phi.html

    static double a1 = 0.254829592;
    static double a2 = -0.284496736;
    static double a3 = 1.421413741;
    static double a4 = -1.453152027;
    static double a5 = 1.061405429;
    static double p = 0.3275911;

    auto sign = x >= 0 ? 1 : -1;
    x = fabs(x) / sqrt(2.0);

    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

    return 0.5 * (1.0 + sign * y);
}

std::vector<double> Construct::DensityEvolutionGA(double snr, size_t nLayers)
{
    static const auto even = [](double x) {
        if (x > 12) {
            return 0.98611 * x - 2.31515;
        }
        if (x > 3.5) {
            return x * (9.0047 * 1e-3 * x + 0.76943) - 0.95068;
        }
        if (x > 1) {
            return x * (0.062883 * x + 0.36784) - 0.16267;
        }
        return x * (0.22024 * x + 0.06448);
    };

    static const auto odd = [](double x) { return x * 2; };

    auto priorMean = std::pow(10, snr / 10) * 2;
    auto means = RunLayersDescent(nLayers, priorMean, even, odd);
    std::transform(means.begin(), means.end(), means.begin(), [](double mean) {
        return 2 * (1 - GaussianCdf(std::sqrt(mean / 2)));
    });

    return means;
}

using Utils::RawPerm;
using Construct::Perm;

std::vector<bool> Construct::BuildFrozenSet_AutomorphismInvariant(
    size_t numLayers, size_t dimension,
    const std::vector<size_t>& blockSizes,
    const std::vector<double>& errorProbs,
    size_t minWeight)
{
    auto length = 1ull << numLayers;
    auto rawPerms = Utils::ListBlockDigitsPermutations(blockSizes);

    std::vector<size_t> lowWeightSymbols;
    std::vector<std::vector<size_t>> classes;
    std::vector classified(length, false);
    std::vector<Perm> perms(rawPerms.size());
    std::transform(rawPerms.begin(), rawPerms.end(), perms.begin(), Perm::MakeDigits);

    for (size_t i = 0; i < length; i++) {
        if (Utils::HammingWeight(i) < minWeight) {
            lowWeightSymbols.push_back(i);
            continue;
        }
        if (classified[i]) {
            continue;
        }

        std::unordered_set<size_t> newClass;
        for (auto& perm : perms) {
            auto j = perm[i];
            newClass.insert(j);
            classified[j] = true;
        }
        classes.emplace_back(newClass.begin(), newClass.end());
    }

    auto frozen = BuildFrozenSet_KnapsackProblemSolution(
        length, dimension + lowWeightSymbols.size(), classes, errorProbs
    );
    for (auto i : lowWeightSymbols) {
        frozen[i] = true;
    }

    return frozen;
}
