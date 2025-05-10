#include <algorithm>
#include <random>
#include <numeric>
#include <vector>
#include <optional>

#include "math/MathUtils.h"
#include "construct/PermSet.h"

using math::Perm;
using math::MatGF2;
using math::VecGF2;

class RandomPermGenerator {
public:
    RandomPermGenerator(size_t nDigits);
    std::vector<size_t> Next();

private:
    size_t mNumDigits;
    std::uniform_int_distribution<size_t> mIndexDistr;
    std::vector<bool> mUsed;
};

RandomPermGenerator::RandomPermGenerator(size_t nDigits)
    : mNumDigits(nDigits)
    , mIndexDistr(0, nDigits - 1)
    , mUsed(nDigits)
{}

std::vector<size_t> RandomPermGenerator::Next()
{
    static std::random_device device;
    static std::mt19937 gen(device());

    std::fill(mUsed.begin(), mUsed.end(), false);
    std::vector<size_t> permutation(mNumDigits);

    for (size_t i = 0; i < mNumDigits; i++)
    {
        size_t next;
        do {
            next = mIndexDistr(gen);
        } while (mUsed[next]);

        mUsed[next] = true;
        permutation[i] = next;
    }

    return permutation;
}

template <typename Container>
static inline size_t HammingDistance(const Container& lhs, const Container& rhs)
{
    size_t count = 0;
    std::for_each(lhs.begin(), lhs.end(), [&, i = 0](auto val) mutable {
        if (val != rhs[i++]) {
            count++;
        }
    });
    return count;
}

std::vector<Perm> construct::BuildDigitsPermSet_Random(size_t count, size_t numLayers, size_t numFrozenLayers, size_t minDist)
{
    RandomPermGenerator generator(numLayers - numFrozenLayers);
    auto id = utils::Iota(numLayers);

    std::vector<std::vector<size_t>> digitPerms = { id };

    std::vector<size_t> perm(numLayers);
    std::iota(perm.begin(), perm.begin() + numFrozenLayers, 0);

    for (size_t i = 0; i < count - 1; i++) {
        auto tail = generator.Next();
        for (size_t j = numFrozenLayers; j < numLayers; j++) {
            perm[j] = tail[j - numFrozenLayers] + numFrozenLayers;
        }
        if (minDist == 0 ||
            std::all_of(digitPerms.begin(), digitPerms.end(),
                [&](const std::vector<size_t>& current) {
                    return HammingDistance(perm, current) >= minDist;
                }
            )) {
            digitPerms.emplace_back(perm);
        }
    }

    std::vector<Perm> result(digitPerms.size());
    std::transform(digitPerms.begin(), digitPerms.end(), result.begin(), Perm::MakeDigits);

    return result;
}

static std::vector<Perm> ListDigitsPerms(size_t numDigits)
{
    auto digitsPerms = utils::ListPermutations(numDigits);
    std::vector<Perm> result(digitsPerms.size());
    std::transform(digitsPerms.begin(), digitsPerms.end(), result.begin(), Perm::MakeDigits);
    return result;
}

static double EstimateErrorProbability(
    const std::vector<double>& errorProbs,
    const std::vector<bool>& frozen)
{
    double successProb = 1;
    for (size_t i = 0; i < errorProbs.size(); i++) {
        if (!frozen[i]) {
            successProb *= (1 - errorProbs[i]);
        }
    }
    return 1 - successProb;
}

static double EstimateErrorProbability_DigitsPerm(
    const Perm& perm,
    const std::vector<double>& errorProbs,
    const std::vector<bool>& frozen)
{
    return EstimateErrorProbability(errorProbs, perm.ApplyDirect(frozen));
}

static codec::PolarSpecification PermuteSpecification_DigitsPerm(
    const codec::PolarSpecification& spec, const Perm& perm)
{
    MatGF2 freezingMatrix = spec.BuildFreezingMatrix();
    MatGF2 permutedFreezingMatrix = MatGF2::Zeros(freezingMatrix.NumRows(), freezingMatrix.NumCols());
    for (size_t i = 0; i < freezingMatrix.NumRows(); i++) {
        for (size_t j = 0; j < freezingMatrix.NumCols(); j++) {
            permutedFreezingMatrix.Set(i, j, freezingMatrix[i][perm[j]]);
        }
    }

    return permutedFreezingMatrix;
}

static double EstimateErrorProbability_Perm(
    const std::vector<double>& errorProbs,
    const MatGF2& permMatrix,
    const MatGF2& freezingMatrix)
{
    codec::PolarSpecification spec = math::TransformFreezingMatrix(freezingMatrix, permMatrix);
    std::vector<bool> frozen(spec.Length);

    for (size_t i = 0; i < spec.Length; i++) {
        frozen[i] = spec.StaticFrozen[i] || spec.Dynamic.Frozen[i];
    }

    return EstimateErrorProbability(errorProbs, frozen);
}

std::vector<Perm> construct::BuildDigitsPermSet_LeastErrorProb(
    size_t count, size_t nLayers, size_t minDist,
    const std::vector<double>& errorProbs,
    const std::vector<bool>& frozen
) {
    auto all = ListDigitsPerms(nLayers);
    std::vector<double> permErrorProbs(all.size());
    std::transform(all.begin(), all.end(), permErrorProbs.begin(), [&](const auto& perm) {
        return EstimateErrorProbability_DigitsPerm(perm, errorProbs, frozen);
        });

    auto indices = utils::SortingPerm(permErrorProbs);
    std::vector<Perm> picked;

    for (auto i : indices) {
        auto& current = all[i];
        auto hasConflict = std::any_of(
            picked.begin(), picked.end(),
            [&](const Perm& perm) {
                return HammingDistance(current.AsVector(), perm.AsVector()) < minDist;
            }
        );

        if (hasConflict) {
            continue;
        }

        picked.push_back(current);

        if (picked.size() == count) {
            return picked;
        }
    }
}

std::vector<Perm> construct::BuildDigitsPermSet_LeastErrorProb(
    size_t count, size_t minDist,
    const std::vector<double>& errorProbs,
    const std::vector<size_t>& blockSizes,
    const codec::PolarSpecification& spec)
{
    auto raw = utils::ListBlockDigitsPermutations(blockSizes);
    std::vector<Perm> all(raw.size());
    std::transform(raw.begin(), raw.end(), all.begin(), Perm::MakeDigits);

    std::vector<double> permErrorProbs(all.size());
    std::transform(
        all.begin(), all.end(), permErrorProbs.begin(),
        [&](const Perm& perm) {
            auto permutedSpec = PermuteSpecification_DigitsPerm(spec, perm);
            double successProb = 1;
            for (size_t i = 0; i < spec.Length; i++) {
                if (!permutedSpec.StaticFrozen[i] && !permutedSpec.Dynamic.Frozen[i]) {
                    successProb *= (1 - errorProbs[i]);
                }
            }
            return 1 - successProb;
        }
    );

    auto indices = utils::SortingPerm(permErrorProbs);
    std::vector<Perm> picked;

    for (auto i : indices) {
        auto& current = all[i];
        auto hasConflict = std::any_of(
            picked.begin(), picked.end(),
            [&](const Perm& perm) {
                return HammingDistance(current.AsVector(), perm.AsVector()) < minDist;
            }
        );

        if (hasConflict) {
            continue;
        }

        picked.push_back(current);

        if (picked.size() == count) {
            return picked;
        }
    }
}

static Perm RandomTrailingDiagPerm(size_t numLayers, size_t blockSize)
{
    auto matrix = MatGF2::Zeros(numLayers, numLayers);
    size_t offset = numLayers - blockSize;

    for (size_t col = 0; col < offset; col++) {
        matrix.Set(col, col, 1);
        /*for (size_t row = col + 1; row < numLayers; row++) {
            matrix.put(row, col, NTL::random_GF2());
        }*/
    }

    MatGF2 trailingBlock;
    bool isLTA;
    do {
        trailingBlock = math::RandomMatGF2(blockSize, blockSize);

        isLTA = true;
        for (size_t row = 0; row < blockSize; row++) {
            for (size_t col = row + 1; col < blockSize; col++) {
                if (trailingBlock[row][col]) {
                    isLTA = false;
                }
            }
        }
    } while (isLTA || trailingBlock.Determinant() == 0);

    for (size_t row = 0; row < blockSize; row++) {
        for (size_t col = 0; col < blockSize; col++) {
            matrix.Set(offset + row, offset + col, trailingBlock[row][col]);
        }
    }

    return Perm::MakeAffine(matrix, VecGF2(numLayers, 0));
}

std::vector<Perm> construct::BuildJointPermSet(size_t count, size_t numLayers, size_t l, size_t h)
{
    auto length = 1ull << numLayers;
    auto blockSize = 1ull << l;

    auto kernel = math::BuildArikanKernel(l);
    auto blockPerms = utils::ListPermutations(blockSize);
    kernel.SelfTranspose();

    static std::random_device device;
    static std::mt19937 gen(device());

    std::uniform_int_distribution<size_t> blockPermDistr(0, blockPerms.size() - 1);
    std::vector<Perm> perms = { Perm::MakeIdentity(length) };
    std::vector<MatGF2> usedMatrices = { MatGF2::Identity(blockSize) };

    size_t d = (l == 0 ? 1 : 4);
    for (size_t i = 0; i < count; i += d) {
        auto tailPerm = RandomTrailingDiagPerm(numLayers, h);

        if (i != 0) {
            std::vector<size_t> indices(length);
            for (size_t b = 0; b < length; b += blockSize) {
                for (size_t j = 0; j < blockSize; j++) {
                    indices[b + j] = tailPerm[b + j];
                }
            }
            perms.push_back(std::move(indices));
        }

        if (l == 0) {
            continue;
        }

        for (size_t p = 0; p < 3; p++) {
            auto& blockPerm = blockPerms[blockPermDistr(gen)];
            auto blockPermMatrix = Perm(blockPerm).AsMatGF2();
            auto matrix = kernel * blockPermMatrix * kernel;
            auto hasConflict = std::any_of(
                usedMatrices.begin(), usedMatrices.end(),
                [&](const MatGF2& usedMatrix) {
                    auto diff = matrix.Invert() * usedMatrix;
                    for (size_t i = 0; i < diff.NumRows(); i++) {
                        for (size_t j = i + 1; j < diff.NumCols(); j++) {
                            if (diff[i][j]) {
                                return false;
                            }
                        }
                    }
                    return true;
                }
            );

            if (l == 0 || !hasConflict) {
                usedMatrices.push_back(std::move(blockPermMatrix));
            }
            else {
                p--;
                continue;
            }

            std::vector<size_t> indices(length);
            for (size_t b = 0; b < length; b += blockSize) {
                for (size_t j = 0; j < blockSize; j++) {
                    indices[b + j] = tailPerm[b + blockPerm[j]];
                }
            }
            perms.push_back(std::move(indices));
        }
    }

    return perms;
}

static VecGF2 GetPolarGeneratorMatrixRow(size_t i, size_t numDigits)
{
    size_t length = 1ull << numDigits;
    VecGF2 vec(length, 0);
    auto digits = utils::Iota(numDigits);

    for (size_t j = 0; j < length; j++) {
        if (std::all_of(digits.begin(), digits.end(), [&](size_t k) {
            auto mask = 1ull << k;
            auto res = (i & mask) != 0 || (j & mask) == 0;
            return res;
        })) {
            vec.Set(j, 1);
        }
    }

    return vec;
}

static MatGF2 BuildPolarSubcodeGeneratorMatrix(const codec::PolarSpecification& spec)
{
    auto numDigits = utils::IntLog2(spec.Length);
    auto mat = MatGF2::Zeros(spec.Dimension, spec.Length);

    std::vector<std::optional<VecGF2>> rows(spec.Length);
    auto GetOrComputeRow = [&](size_t i) -> const VecGF2& {
        if (auto& row = rows[i]) {
            return row.value();
        }
        rows[i] = GetPolarGeneratorMatrixRow(i, numDigits);
        return rows[i].value();
    };

    for (size_t i = 0, r = 0; i < spec.Length; i++) {
        if (spec.StaticFrozen[i] || spec.Dynamic.Frozen[i]) {
            continue;
        }

        VecGF2 row = GetOrComputeRow(i);
        for (size_t dyn : spec.Dynamic.ForwardEquations[i]) {
            row += GetOrComputeRow(dyn);
        }

        for (size_t j = 0; j < spec.Length; j++) {
            mat.Set(r, j, row[j]);
        }

        r++;
    }

    return mat;
}

std::vector<size_t> construct::GetStabBlocksStructure(const std::vector<bool>& frozen)
{
    auto numDigits = utils::IntLog2(frozen.size());

    auto digits = utils::Iota(numDigits);
    auto indices = utils::Iota(frozen.size());

    std::vector<size_t> result;

    for (size_t i = 0; i < numDigits; i++) {
        for (size_t j = numDigits; j-- > i;) {
            digits[i] = j; digits[j] = i;
            auto perm = Perm::MakeDigits(digits);
            digits[i] = i; digits[j] = j;

            if (std::all_of(indices.begin(), indices.end(), [&](size_t k) {
                return frozen[k] == frozen[perm[k]];
            })) {
                result.push_back(j - i + 1);
                i = j;
                break;
            }
        }
    }

    return result;
}

static MatGF2 RandomNonsingularMatGF2(long n)
{
    while (true) {
        auto mat = math::RandomMatGF2(n, n);
        if (mat.Determinant() == 1) {
            return mat;
        }
    }
}

static bool IsSCAbsorbed(const MatGF2& factor)
{
    if (*factor[0].End() > 1) {
        return false;
    }

    for (size_t i = 1; i < factor.NumRows(); i++) {
        if (*factor[i].End() > i) {
            return false;
        }
    }

    return true;
}

std::vector<Perm> construct::BuildBLTAPermSet(
    size_t count, size_t numLayers,
    size_t minDistUTElems, size_t minDistPerm,
    const std::vector<size_t>& blockSizes)
{
    static std::random_device device;
    static std::mt19937 gen(device());

    auto length = 1ull << numLayers;
    auto numBlocks = blockSizes.size();
    auto rawDigitsPerms = utils::ListBlockDigitsPermutations(blockSizes);
    std::uniform_int_distribution<size_t> permIdxDistr(0, rawDigitsPerms.size() - 1);

    size_t numUTElems = std::accumulate(
        blockSizes.begin(), blockSizes.end(), 0,
        [](size_t current, size_t blockSize) {
            return current + blockSize * (blockSize - 1) / 2;
        }
    );

    std::vector<Perm> result = { Perm::MakeIdentity(length) };
    std::vector<utils::RawPerm> usedDigitsPerms = { utils::Iota(numLayers) };
    std::vector<MatGF2> usedFactors = { MatGF2::Identity(numLayers) };
    std::vector<VecGF2> usedUTElems = { VecGF2(numUTElems, 0) };

    while (result.size() < count) {
        auto factor = MatGF2::Identity(numLayers);

        VecGF2 UTElems(numUTElems, 0);
        do {
            UTElems = math::RandomVecGF2(length);
        } while (std::any_of(
            usedUTElems.begin(), usedUTElems.end(),
            [&](const VecGF2& used) {
                return math::HammingDistance(UTElems, used) < minDistUTElems;
            }));

        size_t UTElemIdx = 0;
        size_t rectSize = 0;
        for (auto blockSize : blockSizes) {
            for (size_t i = 0; i < blockSize; i++) {
                for (size_t j = i + 1; j < blockSize; j++) {
                    factor.Set(rectSize + i, rectSize + j, UTElems[UTElemIdx++]);
                }
            }
            rectSize += blockSize;
        }

        while (true) {
            auto& digitsPerm = rawDigitsPerms[permIdxDistr(gen)];
            if (std::any_of(
                usedDigitsPerms.begin(), usedDigitsPerms.end(),
                [&](const utils::RawPerm& used) {
                    return HammingDistance(digitsPerm, used) < minDistPerm;
                }
            )) {
                continue;
            }

            auto permutedFactor = MatGF2::Zeros(numLayers, numLayers);
            for (size_t i = 0; i < numLayers; i++) {
                for (size_t j = 0; j < numLayers; j++) {
                    permutedFactor.Set(i, j, factor[digitsPerm[i]][j]);
                }
            }
            auto invPermutedFactor = permutedFactor.Invert();
            if (std::all_of(
                usedFactors.begin(), usedFactors.end(),
                [&](const MatGF2& usedFactor) {
                    return !IsSCAbsorbed(usedFactor * invPermutedFactor);
                }
            )) {
                factor = std::move(permutedFactor);
                break;
            }
        }

        result.push_back(Perm::MakeAffine(factor, math::RandomVecGF2(numLayers)));
        usedFactors.push_back(std::move(factor));
        usedUTElems.push_back(std::move(UTElems));
    }

    return result;
}
