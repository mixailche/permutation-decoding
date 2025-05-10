#pragma once

#include <cstddef>
#include <numeric>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <array>

namespace utils {
    template <typename T>
    using pair = std::array<T, 2>;

    size_t IntLog2(size_t arg);
    size_t HammingWeight(size_t x);
    bool RandomGF2();

    std::vector<size_t> Iota(size_t n);

    double NegativeMerge(double a, double b);
    double PositiveMerge(double a, double b, bool u);

    using RawPerm = std::vector<size_t>;
    std::vector<RawPerm> ListPermutations(size_t length);

    std::vector<RawPerm> ListBlockDigitsPermutations(const std::vector<size_t>& blockSizes);

    using Comb = uint64_t;
    std::vector<Comb> ListCombinations(size_t length, size_t weight);

    template <typename T, typename Comp = std::less<T>>
    std::vector<size_t> SortingPerm(const std::vector<T>& vec, Comp comp = std::less<T>())
    {
        auto indices = utils::Iota(vec.size());
        std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
            return comp(vec[i], vec[j]);
        });
        return indices;
    }

    template <typename Comp = std::less<double>>
    std::vector<size_t> SolveKnapsackProblem(
        size_t numElements, size_t maxWeight,
        const std::vector<size_t>& weights, const std::vector<double>& costs,
        Comp comp = std::less<double>())
    {
        if (numElements != weights.size() || numElements != costs.size()) {
            throw std::invalid_argument("weights.size() and costs.size() must equal to numElements");
        }

        struct Decision {
            bool Accept;
            double Cost;
        };

        std::vector dp(numElements + 1, std::vector<Decision>(maxWeight + 1, { false, 0 }));

        for (size_t i = 1; i <= numElements; i++) {
            for (size_t totalWeight = 0; totalWeight <= maxWeight; totalWeight++) {
                dp[i][totalWeight] = Decision{ false, dp[i - 1][totalWeight].Cost };
                auto newWeight = weights[i - 1];

                if (newWeight > totalWeight) {
                    continue;
                }

                auto mergedCost = dp[i - 1][totalWeight - newWeight].Cost + costs[i - 1];
                if (comp(mergedCost, dp[i][totalWeight].Cost)) {
                    dp[i][totalWeight] = Decision{ true, mergedCost };
                }
            }
        }

        size_t currentWeight = 0;
        double bestCost = 0;
        for (size_t totalWeight = 0; totalWeight <= maxWeight; totalWeight++) {
            auto cost = dp[numElements][totalWeight].Cost;
            if (comp(cost, bestCost)) {
                bestCost = cost;
                currentWeight = totalWeight;
            }
        }

        std::vector<size_t> result;
        for (size_t i = numElements; i > 0; i--) {
            if (!dp[i][currentWeight].Accept) {
                continue;
            }
            currentWeight -= weights[i - 1];
            result.push_back(i - 1);
        }

        return result;
    }
} // namespace utils
