#include <stdexcept>
#include <random>

#include "Utils.hpp"

size_t utils::IntLog2(size_t arg)
{
    size_t result = -1;

    while (arg > 0) {
        arg /= 2;
        result++;
    };

    return result;
}

size_t utils::HammingWeight(size_t x)
{
    size_t w = 0;
    for (size_t mask = 1; mask != 0; mask <<= 1) {
        if ((x & mask) != 0) {
            w++;
        }
    }
    return w;
}

bool utils::RandomGF2()
{
    static std::mt19937 gen{ std::random_device()() };
    static std::uniform_int_distribution<int> distr(0, 1);

    return distr(gen);
}

std::vector<size_t> utils::Iota(size_t n)
{
    std::vector<size_t> result(n);
    std::iota(result.begin(), result.end(), 0);
    return result;
}

static inline int sgn(double x)
{
    return x > 0 ? 1 : -1;
}

double utils::NegativeMerge(double a, double b)
{
    return sgn(a) * sgn(b) * std::min(std::abs(a), std::abs(b));
}

double utils::PositiveMerge(double a, double b, bool u)
{
    return u == 0 ? a + b : b - a;
}

using utils::RawPerm;

static void ListPermutationsImpl(
    size_t length, RawPerm& current, std::vector<RawPerm>& out,
    std::vector<bool>& used)
{
    if (current.size() == length) {
        return (void)out.push_back(current); // copy
    }

    for (size_t i = 0; i < length; i++) {
        if (used[i]) {
            continue;
        }

        current.push_back(i);
        used[i] = true;
        ListPermutationsImpl(length, current, out, used);
        current.pop_back();
        used[i] = false;
    }
}

std::vector<RawPerm> utils::ListPermutations(size_t length)
{
    std::vector<RawPerm> out;
    std::vector<bool> used(length, false);
    RawPerm current = {};
    ListPermutationsImpl(length, current, out, used);

    return out;
}

static void ListBlockDigitsPermutationsImpl(
    const std::vector<size_t>& blockSizes,
    size_t currentBlockNumber, size_t offset,
    RawPerm& currentPerm, std::vector<RawPerm>& out)
{
    if (currentBlockNumber == blockSizes.size()) {
        return (void)out.push_back(currentPerm); // copy
    }

    auto blockSize = blockSizes[currentBlockNumber];
    auto tails = utils::ListPermutations(blockSize);

    for (auto& tail : tails) {
        for (size_t i = 0; i < blockSize; i++) {
            currentPerm[offset + i] = offset + tail[i];
        }
        ListBlockDigitsPermutationsImpl(
            blockSizes, currentBlockNumber + 1, offset + blockSize, currentPerm, out
        );
    }
}

std::vector<RawPerm> utils::ListBlockDigitsPermutations(const std::vector<size_t>& blockSizes)
{
    RawPerm currentPerm(std::accumulate(blockSizes.begin(), blockSizes.end(), 0ull));
    std::vector<RawPerm> out;
    ListBlockDigitsPermutationsImpl(blockSizes, 0, 0, currentPerm, out);
    return out;
}

using utils::Comb;

static void ListCombinationsImpl(
    size_t length, size_t weight, size_t pos, size_t currentWeight,
    Comb currentComb, std::vector<Comb>& out)
{
    if (pos == length) {
        return out.push_back(currentComb);
    }

    if (weight - currentWeight < length - pos) {
        ListCombinationsImpl(length, weight, pos + 1, currentWeight, currentComb, out);
    }

    if (currentWeight < weight) {
        ListCombinationsImpl(length, weight, pos + 1, currentWeight + 1,
            currentComb | (1ull << pos), out);
    }
}

std::vector<Comb> utils::ListCombinations(size_t length, size_t weight)
{
    std::vector<Comb> out;
    ListCombinationsImpl(length, weight, 0, 0, 0, out);
    return out;
}
