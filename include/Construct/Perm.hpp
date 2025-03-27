#pragma once

#include <vector>

#include "Specification.h"

#include "NTLHelp.h"
#include "NTL/GF2E.h"
#include "NTL/mat_GF2.h"
#include "NTL/vec_GF2.h"

namespace Construct {
    class Perm {
    public:
        Perm() = default;
        Perm(std::vector<size_t> newIndices);

        static Perm MakeIdentity(size_t length);
        static Perm MakeDigits(const std::vector<size_t>& digitsPerm);
        static Perm MakeAffine(NTL::mat_GF2 factor, NTL::vec_GF2 shift);

        size_t Length() const;

        size_t& operator[](size_t number);
        const size_t& operator[](size_t number) const;

        const std::vector<size_t>& AsVector() const;
        const NTL::mat_GF2 AsMatGF2() const;

        template <typename T>
        std::vector<T> ApplyDirect(const std::vector<T>& vec) const;
    
        template <typename T>
        std::vector<T> ApplyReversed(const std::vector<T>& vec) const;

    private:
        std::vector<size_t> mNewIndices;
    };

    Codec::PolarSubcodeSpecification PermuteSpecification(
        const Codec::PolarSubcodeSpecification& spec,
        const Perm& perm
    );

    std::vector<size_t> GetStabBlocksStructure(const std::vector<bool>& frozen);

    std::vector<size_t> GetStabBlocksStructure(
        const std::vector<bool>& frozen,
        const std::vector<double>& errorProbs, double alpha
    );

    std::vector<Perm> BuildDigitsPermSet_Random(
        size_t count, size_t numLayers, size_t numFrozenLayers, size_t minDist
    );

    std::vector<Perm> BuildDigitsPermSet_LeastErrorProb(
        size_t count, size_t numLayers, size_t minDist,
        const std::vector<double>& errorProbs,
        const std::vector<bool>& frozen
    );

    std::vector<Perm> BuildDigitsPermSet_LeastErrorProb(
        size_t count, size_t minDist,
        const std::vector<double>& errorProbs,
        const std::vector<size_t>& blockSizes,
        const Codec::PolarSubcodeSpecification& spec
    );

    std::vector<Perm> BuildBlockPermSet_LeastErrorProb(
        size_t count, size_t numLayers,
        size_t l, size_t h
    );

    std::vector<Perm> BuildAffineGF2EPermSet_LeastErrorProb(
        size_t count, size_t minDistFactor, size_t minDistShift,
        const std::vector<double>& errorProbs,
        const Codec::PolarSubcodeSpecification& spec
    );

    std::vector<Perm> BuildBLTAPermSet(
        size_t count, size_t numLayers,
        size_t minDistUTElems, size_t minDistPerm,
        const std::vector<size_t>& blockSizes
    );

    std::vector<Perm> BuildBLTAPermSet(
        size_t count, size_t minDistUTElems, size_t minDistPerm,
        const std::vector<size_t>& blockSizes,
        const std::vector<double>& errorProbs,
        const Codec::PolarSubcodeSpecification& spec
    );

    std::vector<Perm> BuildRMPolarPermSet(
        size_t count, size_t minDistUTElems, size_t minDistPerm,
        const Codec::PolarCodeSpecification& outerPolarSpec
    );
} // namespace construct

using Construct::Perm;

template <typename T>
std::vector<T> Perm::ApplyDirect(const std::vector<T>& vec) const
{
    std::vector<T> result(vec.size());

    for (size_t i = 0; i < vec.size(); i++) {
        result[i] = vec[this->operator[](i)];
    }

    return result;
}

template <typename T>
std::vector<T> Perm::ApplyReversed(const std::vector<T>& vec) const
{
    std::vector<T> result(vec.size());

    for (size_t i = 0; i < vec.size(); i++) {
        result[this->operator[](i)] = vec[i];
    }

    return result;
}
