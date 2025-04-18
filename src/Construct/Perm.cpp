#include <algorithm>
#include <random>
#include <numeric>
#include <vector>
#include <optional>

#include "NTLHelp.h"
#include "Utils.hpp"
#include "Specification.h"
#include "Construct/Perm.hpp"
#include "Construct/PolarSubcode.h"

using Construct::Perm;

Perm::Perm(std::vector<size_t> newIndices)
    : mNewIndices(std::move(newIndices))
{}

Perm Perm::MakeIdentity(size_t length)
{
    return Utils::Iota(length);
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

Perm Perm::MakeDigits(const std::vector<size_t>& digitsPerm)
{
    std::vector<size_t> newIndices(1ull << digitsPerm.size());
    for (size_t i = 0; i < newIndices.size(); i++) {
        newIndices[i] = ApplyDigitsPerm(digitsPerm, i);
    }
    return std::move(newIndices);
}

Perm Perm::MakeAffine(NTL::mat_GF2 factor, NTL::vec_GF2 shift)
{
    size_t length = factor.NumRows();
    size_t cardinality = 1ull << length;
    std::vector<size_t> newIndices(cardinality);

    for (size_t i = 0; i < cardinality; i++) {
        auto digitsVector = Math::IndexToGF2Vec(i, length);
        newIndices[i] = Math::GF2VecToIndex(factor * digitsVector + shift);
    }

    return newIndices;
}

size_t Perm::Length() const
{
    return mNewIndices.size();
}

size_t& Perm::operator[](size_t number)
{
    return mNewIndices[number];
}

const size_t& Perm::operator[](size_t number) const
{
    return mNewIndices[number];
}

const std::vector<size_t>& Perm::AsVector() const
{
    return mNewIndices;
}

const NTL::mat_GF2 Perm::AsMatGF2() const
{
    NTL::mat_GF2 permMatrix(NTL::INIT_SIZE, Length(), Length());
    for (size_t i = 0; i < mNewIndices.size(); i++) {
        permMatrix.put(i, mNewIndices[i], 1);
    }
    return permMatrix;
}

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

template <typename Container> static inline size_t
HammingDistance(const Container& lhs, const Container& rhs)
{
    size_t count = 0;
    std::for_each(lhs.begin(), lhs.end(), [&, i = 0](auto val) mutable {
        if (val != rhs[i++]) {
            count++;
        }
    });
    return count;
}

std::vector<Perm> Construct::BuildDigitsPermSet_Random(
    size_t count, size_t numLayers, size_t numFrozenLayers, size_t minDist)
{
    RandomPermGenerator generator(numLayers - numFrozenLayers);
    auto id = Utils::Iota(numLayers);

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
    auto digitsPerms = Utils::ListPermutations(numDigits);
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

static Codec::PolarSubcodeSpecification PermuteFreezingMatrix(
    const NTL::mat_GF2& freezingMatrix, const Perm& perm)
{
    return Math::TransformFreezingMatrix_Perm(freezingMatrix, perm.AsMatGF2());
}

Codec::PolarSubcodeSpecification Construct::PermuteSpecification(
    const Codec::PolarSubcodeSpecification& spec,
    const Perm& perm)
{
    return PermuteFreezingMatrix(BuildFreezingMatrix(spec), perm);
}

static Codec::PolarSubcodeSpecification PermuteSpecification_DigitsPerm(
    const Codec::PolarSubcodeSpecification& spec,
    const Perm& perm)
{
    NTL::mat_GF2 freezingMatrix = Construct::BuildFreezingMatrix(spec);
    NTL::mat_GF2 permutedFreezingMatrix(NTL::INIT_SIZE, freezingMatrix.NumRows(), freezingMatrix.NumCols());
    for (size_t i = 0; i < freezingMatrix.NumRows(); i++) {
        for (size_t j = 0; j < freezingMatrix.NumCols(); j++) {
            permutedFreezingMatrix[i][j] = freezingMatrix[i][perm[j]];
        }
    }

    Codec::PolarSubcodeSpecification permutedSpec(spec.Length);
    permutedSpec.Dimension = spec.Dimension;

    Math::RunDoubleGaussianElimination(permutedFreezingMatrix);
    auto equations = Construct::BuildDynamicFreezingEquations(permutedFreezingMatrix, permutedSpec.StaticFrozen);
    for (size_t i = 0; i < spec.Length; i++) {
        for (auto j : equations[i]) {
            permutedSpec.Dynamic.Frozen[i] = true;
            permutedSpec.Dynamic.ForwardEquations[j].push_back(i);
        }
    }

    return permutedSpec;
}

static double EstimateErrorProbability_Perm(
    const NTL::mat_GF2& permMatrix,
    const std::vector<double>& errorProbs,
    const NTL::mat_GF2& freezingMatrix)
{
    Codec::PolarSubcodeSpecification spec = Math::TransformFreezingMatrix_Perm(freezingMatrix, permMatrix);
    std::vector<bool> frozen(spec.Length);

    for (size_t i = 0; i < spec.Length; i++) {
        frozen[i] = spec.StaticFrozen[i] || spec.Dynamic.Frozen[i];
    }
    
    return EstimateErrorProbability(errorProbs, frozen);
}

std::vector<Perm> Construct::BuildDigitsPermSet_LeastErrorProb(
    size_t count, size_t nLayers, size_t minDist,
    const std::vector<double>& errorProbs,
    const std::vector<bool>& frozen
) {
    auto all = ListDigitsPerms(nLayers);
    std::vector<double> permErrorProbs(all.size());
    std::transform(all.begin(), all.end(), permErrorProbs.begin(), [&](const auto& perm) {
        return EstimateErrorProbability_DigitsPerm(perm, errorProbs, frozen);
    });

    auto indices = Utils::SortingPerm(permErrorProbs);
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

std::vector<Perm> Construct::BuildDigitsPermSet_LeastErrorProb(
    size_t count, size_t minDist,
    const std::vector<double>& errorProbs,
    const std::vector<size_t>& blockSizes,
    const Codec::PolarSubcodeSpecification& spec)
{
    auto raw = Utils::ListBlockDigitsPermutations(blockSizes);
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

    auto indices = Utils::SortingPerm(permErrorProbs);
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
    NTL::mat_GF2 matrix(NTL::INIT_SIZE, numLayers, numLayers);
    size_t offset = numLayers - blockSize;

    for (size_t col = 0; col < offset; col++) {
        matrix.put(col, col, 1);
        /*for (size_t row = col + 1; row < numLayers; row++) {
            matrix.put(row, col, NTL::random_GF2());
        }*/
    }

    NTL::mat_GF2 trailingBlock;
    bool isLTA;
    do {
        trailingBlock = NTL::random_mat_GF2(blockSize, blockSize);

        isLTA = true;
        for (size_t row = 0; row < blockSize; row++) {
            for (size_t col = row + 1; col < blockSize; col++) {
                if (NTL::IsOne(trailingBlock[row][col])) {
                    isLTA = false;
                }
            }
        }
    } while (isLTA || NTL::IsZero(NTL::determinant(trailingBlock)));

    for (size_t row = 0; row < blockSize; row++) {
        for (size_t col = 0; col < blockSize; col++) {
            matrix[offset + row][offset + col] = trailingBlock[row][col];
        }
    }

    return Perm::MakeAffine(matrix, NTL::vec_GF2(NTL::INIT_SIZE, numLayers));
}

std::vector<Perm> Construct::BuildJointPermSet(
    size_t count, size_t numLayers,
    size_t l, size_t h)
{
    auto length = 1ull << numLayers;
    auto blockSize = 1ull << l;

    auto kernel = NTL::transpose(Math::BuildArikanKernel(l));
    auto blockPerms = Utils::ListPermutations(blockSize);

    static std::random_device device;
    static std::mt19937 gen(device());

    std::uniform_int_distribution<size_t> blockPermDistr(0, blockPerms.size() - 1);
    std::vector<Perm> perms = { Perm::MakeIdentity(length) };
    std::vector<NTL::mat_GF2> usedMatrices = { NTL::ident_mat_GF2(blockSize) };

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
                [&](const NTL::mat_GF2& usedMatrix) {
                    auto diff = NTL::inv(matrix) * usedMatrix;
                    for (size_t i = 0; i < diff.NumRows(); i++) {
                        for (size_t j = i + 1; j < diff.NumCols(); j++) {
                            if (NTL::IsOne(diff[i][j])) {
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

static NTL::vec_GF2 GetPolarGeneratorMatrixRow(size_t i, size_t numDigits)
{
    size_t length = 1ull << numDigits;
    NTL::vec_GF2 vec(NTL::INIT_SIZE, length);
    auto digits = Utils::Iota(numDigits);

    for (size_t j = 0; j < length; j++) {
        if (std::all_of(digits.begin(), digits.end(), [&](size_t k) {
                auto mask = 1ull << k;
                auto res = (i & mask) != 0 || (j & mask) == 0;
                return res;
            })) {
            vec[j] = 1;
        }
    }

    return vec;
}

static NTL::mat_GF2 BuildPolarSubcodeGeneratorMatrix(const Codec::PolarSubcodeSpecification& spec)
{
    auto numDigits = Utils::IntLog2(spec.Length);
    NTL::mat_GF2 mat(NTL::INIT_SIZE, spec.Dimension, spec.Length);

    std::vector<std::optional<NTL::vec_GF2>> rows(spec.Length);
    auto GetOrComputeRow = [&](size_t i) -> const NTL::vec_GF2& {
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
        
        NTL::vec_GF2 row = GetOrComputeRow(i);
        for (size_t dyn : spec.Dynamic.ForwardEquations[i]) {
            row += GetOrComputeRow(dyn);
        }
        
        for (size_t j = 0; j < spec.Length; j++) {
            mat.put(r, j, row[j]);
        }

        r++;
    }

    return mat;
}

std::vector<Perm> Construct::BuildAffineGF2EPermSet_LeastErrorProb(
    size_t count, size_t minDistFactor, size_t minDistShift,
    const std::vector<double>& errorProbs,
    const Codec::PolarSubcodeSpecification& spec)
{
    static auto x = NTL::GF2X(NTL::INIT_MONO, 1, 1);
    static auto primitiveElem = NTL::conv<NTL::GF2E>(x);

    if (count == 0) {
        throw std::invalid_argument("Parameter `count` must be positive");
    }

    size_t deg = NTL::GF2E::degree();
    size_t cardinality = 1ull << deg;
    auto zeroShift = NTL::vec_GF2(NTL::INIT_SIZE, deg);

    //auto freezingMatrix = BuildFreezingMatrix(spec);
    //std::vector<std::pair<Perm, double>> factorPermErrorProbs;

    //for (size_t factorDigits = 1; factorDigits < cardinality; factorDigits++) {
    //    std::cout << "elem = " << factorDigits << std::endl;
    //    NTL::mat_GF2 factor(NTL::INIT_SIZE, deg, deg);
    //    NTL::GF2E factorElem = Math::IndexToGF2E(factorDigits, deg);
    //    
    //    for (size_t i = 0; i < deg; i++) {
    //        auto row = NTL::conv<NTL::GF2X>(factorElem);
    //        factorElem *= primitiveElem;
    //        for (size_t j = 0; j < deg; j++) {
    //            factor[i][j] = NTL::coeff(row, j);
    //        }
    //    }

    //    auto perm = Perm::MakeAffine(factor, zeroShift);
    //    auto prob = EstimateErrorProbability_Perm(perm.AsMatGF2(), errorProbs, freezingMatrix);
    //    factorPermErrorProbs.emplace_back(std::move(perm), prob);
    //}

    //std::sort(factorPermErrorProbs.begin(), factorPermErrorProbs.end(),
    //    [&](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });

    //std::vector<Perm> result(count);
    //std::transform(
    //    factorPermErrorProbs.begin(), factorPermErrorProbs.begin() + count, result.begin(),
    //    [](std::pair<Perm, double>& pair) { return std::move(pair.first); }
    //);
    //return result;

    auto genMatrix = BuildPolarSubcodeGeneratorMatrix(spec);
    auto invGenMatrix = Math::InvRectMatrix(genMatrix);

    std::vector<NTL::vec_GF2> usedFactors;
    std::vector<NTL::vec_GF2> usedShifts;
    auto HasConflict = [](
        const std::vector<NTL::vec_GF2>& used,
        const NTL::vec_GF2& vec, size_t minDist)
    {
        return std::any_of(used.begin(), used.end(), [&](const NTL::vec_GF2& usedVec) {
            return HammingDistance(vec, usedVec) < minDist;
        });
    };

    std::vector<Perm> result;

    for (size_t factorDigits = 1; factorDigits < cardinality; factorDigits++) {
        NTL::GF2E factorElem = Math::IndexToGF2E(factorDigits, deg);
        NTL::vec_GF2 factorVec = NTL::to_vec_GF2(NTL::conv<NTL::GF2X>(factorElem));
        if (HasConflict(usedFactors, factorVec, minDistFactor)) {
            continue;
        }

        NTL::mat_GF2 factor(NTL::INIT_SIZE, deg, deg);
        
        for (size_t i = 0; i < deg; i++) {
            auto row = NTL::conv<NTL::GF2X>(factorElem);
            for (size_t j = 0; j < deg; j++) {
                factor[i][j] = NTL::coeff(row, j);
            }
            factorElem *= primitiveElem;
        }

        auto dummyPerm = Perm::MakeAffine(factor, zeroShift);
        NTL::mat_GF2 permutedGenMatrix(NTL::INIT_SIZE, spec.Dimension, spec.Length);
        for (size_t i = 0; i < spec.Dimension; i++) {
            for (size_t j = 0; j < spec.Length; j++) {
                permutedGenMatrix[i][j] = genMatrix[i][dummyPerm[j]];
            }
        }

        if (NTL::IsZero(NTL::determinant(permutedGenMatrix * invGenMatrix))) {
            continue;
        }

        std::cout << "elem = " << factorDigits << "\n";

        NTL::vec_GF2 shiftVec;
        for (size_t shiftDigits = 0; shiftDigits < cardinality; shiftDigits++) {
            shiftVec = Math::IndexToGF2Vec(shiftDigits, deg);
            if (!HasConflict(usedShifts, shiftVec, minDistShift)) {
                break;
            }
        }

        result.push_back(Perm::MakeAffine(factor, zeroShift));
        usedFactors.push_back(std::move(factorVec));
        usedShifts.push_back(std::move(shiftVec));

        if (result.size() == count) {
            return result;
        }
    }

    return result;
}

std::vector<size_t> Construct::GetStabBlocksStructure(const std::vector<bool>& frozen)
{
    auto numDigits = Utils::IntLog2(frozen.size());

    auto digits = Utils::Iota(numDigits);
    auto indices = Utils::Iota(frozen.size());

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

std::vector<size_t> Construct::GetStabBlocksStructure(
    const std::vector<bool>& frozen,
    const std::vector<double>& errorProbs, double alpha)
{
    auto numDigits = Utils::IntLog2(frozen.size());
    
    auto digits = Utils::Iota(numDigits);
    auto indices = Utils::Iota(frozen.size());

    std::vector<double> logs(frozen.size());
    std::transform(errorProbs.begin(), errorProbs.end(), logs.begin(), [&](double p) {
        return std::log(p);
    });

    std::vector<size_t> result;

    for (size_t i = 0; i < numDigits; i++) {
        for (size_t j = numDigits; j-- > i;) {
            digits[i] = j; digits[j] = i;
            auto perm = Perm::MakeDigits(digits);
            digits[i] = i; digits[j] = j;

            if (std::all_of(indices.begin(), indices.end(), [&](size_t k) {
                return frozen[k] == frozen[perm[k]]
                    || std::abs(errorProbs[k] - errorProbs[perm[k]]) < alpha;
            })) {
                result.push_back(j - i + 1);
                i = j;
                break;
            }
        }
    }

    return result;
}

static NTL::mat_GF2 RandomNonsingularMatGF2(long n)
{
    while (true) {
        auto mat = NTL::random_mat_GF2(n, n);
        if (NTL::IsOne(NTL::determinant(mat))) {
            return mat;
        }
    }
}

static bool IsSCAbsorbed(const NTL::mat_GF2& factor)
{
    if (*Math::HighestOnePos(factor[0]) > 1) {
        return false;
    }

    for (size_t i = 1; i < factor.NumRows(); i++) {
        if (*Math::HighestOnePos(factor[i]) > i) {
            return false;
        }
    }

    return true;
}

std::vector<Perm> Construct::BuildBLTAPermSet(
    size_t count, size_t numLayers,
    size_t minDistUTElems, size_t minDistPerm,
    const std::vector<size_t>& blockSizes)
{
    static std::random_device device;
    static std::mt19937 gen(device());

    auto length = 1ull << numLayers;
    auto numBlocks = blockSizes.size();
    auto rawDigitsPerms = Utils::ListBlockDigitsPermutations(blockSizes);
    std::uniform_int_distribution<size_t> permIdxDistr(0, rawDigitsPerms.size() - 1);

    size_t numUTElems = std::accumulate(
        blockSizes.begin(), blockSizes.end(), 0,
        [](size_t current, size_t blockSize) {
            return current + blockSize * (blockSize - 1) / 2;
        }
    );

    std::vector<Perm> result = { Perm::MakeIdentity(length) };
    std::vector<Utils::RawPerm> usedDigitsPerms = { Utils::Iota(numLayers) };
    std::vector<NTL::mat_GF2> usedFactors = { NTL::ident_mat_GF2(numLayers) };
    std::vector<NTL::vec_GF2> usedUTElems = { NTL::vec_GF2(NTL::INIT_SIZE, numUTElems) };

    while (result.size() < count) {
        auto factor = NTL::ident_mat_GF2(numLayers);

        NTL::vec_GF2 UTElems(NTL::INIT_SIZE, numUTElems);
        do {
            UTElems = NTL::random_vec_GF2(length);
        } while (std::any_of(
            usedUTElems.begin(), usedUTElems.end(),
            [&](const NTL::vec_GF2& used) {
                return HammingDistance(UTElems, used) < minDistUTElems;
            }));
        
        size_t UTElemIdx = 0;
        size_t rectSize = 0;
        for (auto blockSize : blockSizes) {
            for (size_t i = 0; i < blockSize; i++) {
                for (size_t j = i + 1; j < blockSize; j++) {
                    factor[rectSize + i][rectSize + j] = UTElems[UTElemIdx++];
                }
            }
            rectSize += blockSize;
        }

        while (true) {
            auto& digitsPerm = rawDigitsPerms[permIdxDistr(gen)];
            if (std::any_of(
                usedDigitsPerms.begin(), usedDigitsPerms.end(),
                [&](const Utils::RawPerm& used) {
                    return HammingDistance(digitsPerm, used) < minDistPerm;
                }
            )) {
                continue;
            }

            NTL::mat_GF2 permutedFactor(NTL::INIT_SIZE, numLayers, numLayers);
            for (size_t i = 0; i < numLayers; i++) {
                for (size_t j = 0; j < numLayers; j++) {
                    permutedFactor[i][j] = factor[digitsPerm[i]][j];
                }
            }
            auto invPermutedFactor = NTL::inv(permutedFactor);
            if (std::all_of(
                usedFactors.begin(), usedFactors.end(),
                [&](const NTL::mat_GF2& usedFactor) {
                    return !IsSCAbsorbed(usedFactor * invPermutedFactor);
                }
            )) {
                factor = std::move(permutedFactor);
                break;
            }
        }

        result.push_back(Perm::MakeAffine(factor, NTL::vec_GF2(NTL::INIT_SIZE, numLayers)));
        usedFactors.push_back(std::move(factor));
        usedUTElems.push_back(std::move(UTElems));
    }

    return result;
}

std::vector<Perm> Construct::BuildBLTAPermSet(
    size_t count, size_t minDistUTElems, size_t minDistPerm,
    const std::vector<size_t>& blockSizes,
    const std::vector<double>& errorProbs,
    const Codec::PolarSubcodeSpecification& spec)
{
    static std::random_device device;
    static std::mt19937 gen(device());

    auto freezingMatrix = BuildFreezingMatrix(spec);
    auto numLayers = Utils::IntLog2(spec.Length);
    auto numBlocks = blockSizes.size();
    auto rawDigitsPerms = Utils::ListBlockDigitsPermutations(blockSizes);
    std::uniform_int_distribution<size_t> permIdxDistr(0, rawDigitsPerms.size() - 1);

    /*auto genMatrix = BuildPolarSubcodeGeneratorMatrix(spec);
    auto invGenMatrix = Math::InvRectMatrix(genMatrix);*/

    std::vector rectSizes = { 0 };
    for (size_t i = 1; i < numBlocks; i++) {
        rectSizes.push_back(rectSizes[i - 1] + blockSizes[i - 1]);
    }

    size_t numUTElems = std::accumulate(
        blockSizes.begin(), blockSizes.end(), 0,
        [](size_t current, size_t blockSize) {
            return current + blockSize * (blockSize - 1) / 2;
        }
    );

    struct Entry {
        NTL::vec_GF2 UpperTriangElems;
        NTL::mat_GF2 Factor;
        Utils::RawPerm RawDigitsPerm;
        Perm Perm;
        double ErrorProb;
    };

    std::vector<Entry> entries = { {
        NTL::vec_GF2(NTL::INIT_SIZE, numUTElems),
        NTL::ident_mat_GF2(numLayers),
        Utils::Iota(numLayers),
        Perm::MakeIdentity(spec.Length),
        0 // identity permutation must be included
    } };

    while (entries.size() < 2000) {
        auto factor = NTL::ident_mat_GF2(numLayers);

        auto upperTriangElems = NTL::random_vec_GF2(numUTElems);
        auto rawDigitsPerm = rawDigitsPerms[permIdxDistr(gen)];
        auto digitsPerm = Perm::MakeDigits(rawDigitsPerm);
        
        size_t UTElemIdx = 0;
        size_t rectSize = 0;
        for (auto blockSize : blockSizes) {
            for (size_t i = 0; i < blockSize; i++) {
                for (size_t j = i + 1; j < blockSize; j++) {
                    factor[rectSize + i][rectSize + j] = upperTriangElems[UTElemIdx++];
                }
            }
            rectSize += blockSize;
        }

        NTL::mat_GF2 permutedFactor(NTL::INIT_SIZE, numLayers, numLayers);
        for (size_t i = 0; i < numLayers; i++) {
            for (size_t j = 0; j < numLayers; j++) {
                permutedFactor[i][j] = factor[rawDigitsPerm[i]][j];
            }
        }

        auto perm = Perm::MakeAffine(permutedFactor, NTL::vec_GF2(NTL::INIT_SIZE, numLayers));
        auto errorProb = EstimateErrorProbability_Perm(perm.AsMatGF2(), errorProbs, freezingMatrix);
        
        entries.push_back({
            std::move(upperTriangElems),
            std::move(permutedFactor),
            std::move(rawDigitsPerm),
            std::move(perm),
            errorProb
        });
    }

    std::sort(entries.begin(), entries.end(), [&](const Entry& lhs, const Entry& rhs) {
        return lhs.ErrorProb < rhs.ErrorProb;
    });
    
    std::vector<size_t> indices;
    for (size_t i = 0; i < entries.size() && indices.size() < count; i++) {
        auto& current = entries[i];
        auto invFactor = NTL::inv(current.Factor);
        if (std::all_of(indices.begin(), indices.end(), [&](size_t j) {
            auto& other = entries[j];
            return !IsSCAbsorbed(other.Factor * invFactor)
                && HammingDistance(current.UpperTriangElems, other.UpperTriangElems) >= minDistUTElems
                && HammingDistance(current.RawDigitsPerm, other.RawDigitsPerm) >= minDistPerm;
            })) {
            indices.push_back(i);
        }
    }

    std::vector<Perm> result(indices.size());
    std::transform(indices.begin(), indices.end(), result.begin(),
        [&](size_t i) { return entries[i].Perm; });

    return result;
}

std::vector<Perm> Construct::BuildRMPolarPermSet(
    size_t count, size_t minDistUTElems, size_t minDistPerm,
    const Codec::PolarCodeSpecification& outerPolarSpec)
{
    auto numLayers = Utils::IntLog2(outerPolarSpec.Length);
    auto blockSizes = GetStabBlocksStructure(outerPolarSpec.Frozen);
    
    auto permsEven = BuildBLTAPermSet(count, numLayers, minDistPerm, minDistPerm, blockSizes);
    auto permsOdd = BuildBLTAPermSet(count, numLayers, minDistPerm, minDistPerm, blockSizes);
    std::vector<Perm> perms;

    for (size_t permIdx = 0; permIdx < count; permIdx++) {
        auto& newIndicesEven = permsEven[permIdx].AsVector();
        auto& newIndicesOdd = permsOdd[permIdx].AsVector();

        std::vector<size_t> newIndices(1ull << (numLayers + 1));
        for (size_t i = 0; i < outerPolarSpec.Length; i++) {
            newIndices[2 * i] = 2 * newIndicesEven[i];
            newIndices[2 * i + 1] = 2 * newIndicesOdd[i] + 1;
        }

        perms.emplace_back(std::move(newIndices));
    }

    return perms;
}
