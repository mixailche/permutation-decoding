#include <bitset>

#include "Linalg/BinVector.h"

using Linalg::BinVector;

BinVector::BinVector(size_t dimension)
    : mDimension(dimension)
    , mNBlocks((dimension + BlockSize - 1) / BlockSize)
    , mBlocks(new block_t[mNBlocks])
{}

BinVector::BinVector(size_t dim, bool value)
    : BinVector(dim)
{
    auto block = value ? Ones : 0;

    for (size_t i = 0; i < mNBlocks - 1; i++) {
        mBlocks[i] = block;
    }

    auto nTrailingBits = dim % BlockSize;
    mBlocks[mNBlocks - 1] = nTrailingBits == 0 ? block
        : block & ((1ull << nTrailingBits) - 1);
}

BinVector::BinVector(const std::vector<bool>& symbols)
    : BinVector(symbols.size(), 0)
{
    block_t mask = 1;
    mBlocks[0] = 0;

    for (size_t pos = 0, blockIdx = 0; pos < mDimension; pos++) {
        if (pos > 0 && pos % BlockSize == 0) {
            mask = 1;
            mBlocks[++blockIdx] = 0;
        }

        if (symbols[pos]) {
            mBlocks[blockIdx] |= mask;
        }

        mask <<= 1;
    }
}

BinVector::BinVector(const BinVector& other)
    : BinVector(other.mDimension)
{
    std::copy_n(mBlocks.get(), mNBlocks, other.mBlocks.get());
}

BinVector::BinVector(BinVector&& other) noexcept
    : mNBlocks(other.mNBlocks)
    , mDimension(other.mDimension)
    , mBlocks(other.mBlocks.release())
{}

BinVector& BinVector::operator=(const BinVector& other)
{
    mNBlocks = other.mNBlocks;
    mDimension = other.mDimension;
    mBlocks.reset(new block_t[mNBlocks]);
    std::copy_n(mBlocks.get(), mNBlocks, other.mBlocks.get());
    return *this;
}

BinVector& BinVector::operator=(BinVector&& other) noexcept
{
    mNBlocks = other.mNBlocks;
    mDimension = other.mDimension;
    mBlocks = std::move(other.mBlocks);
    return *this;
}

size_t BinVector::Dimension() const
{
    return mDimension;
}

std::optional<size_t> BinVector::end() const
{
    size_t i = mNBlocks - 1;
    while (i > 0 && mBlocks[i] == 0) {
        i--;
    }

    if (i == 0 && mBlocks[0] == 0)
        return std::nullopt;

    auto block = mBlocks[i];
    size_t highestOnePos = 0;

    while (block > 0) {
        block /= 2;
        highestOnePos++;
    }

    return i * BlockSize + highestOnePos - 1;
}

bool BinVector::operator[](size_t pos) const
{
    auto offset = pos / BlockSize;
    auto index = pos % BlockSize;
    return (mBlocks[offset] & (1ull << index)) != 0;
}

void BinVector::Set(size_t pos, bool value)
{
    auto offset = pos / BlockSize;
    auto index = pos % BlockSize;

    if (value) {
        mBlocks[offset] |= 1ull << index;
    }
    else {
        mBlocks[offset] &= Ones ^ (1ull << index);
    }
}

bool BinVector::operator*(const BinVector& other) const
{
    auto max_dim = std::max(mDimension, other.mDimension);
    BinVector result(max_dim);

    size_t count = 0;

    for (size_t i = 0; i < result.mNBlocks; i++) {
        auto thisBlock = mBlocks[i];
        auto otherBlock = other.mBlocks[i];
        count += std::bitset<BlockSize>(thisBlock & otherBlock).count();
    }

    return count % 2;
}

BinVector BinVector::operator+(const BinVector& other) const
{
    auto maxDimension = std::max(mDimension, other.mDimension);
    BinVector result(maxDimension);

    for (size_t i = 0; i < result.mNBlocks; i++) {
        result.mBlocks[i] = mBlocks[i] ^ other.mBlocks[i];
    }

    return result;
}

BinVector& BinVector::operator+=(const BinVector& other)
{
    for (size_t i = 0; i < other.mNBlocks; i++) {
        mBlocks[i] ^= other.mBlocks[i];
    }

    return *this;
}
