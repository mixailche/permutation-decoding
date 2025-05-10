#include <bitset>

#include "VecGF2.h"

using math::VecGF2;

VecGF2::VecGF2(size_t dimension)
    : mDimension(dimension)
    , mNumBlocks((dimension + BLOCK_SIZE - 1) / BLOCK_SIZE)
    , mBlocks(new block_t[mNumBlocks])
{}

VecGF2::VecGF2()
    : VecGF2(0)
{}

VecGF2::VecGF2(size_t dim, bool value)
    : VecGF2(dim)
{
    if (dim == 0) {
        return;
    }

    auto block = value ? ONES : 0;

    for (size_t i = 0; i < mNumBlocks - 1; i++) {
        mBlocks[i] = block;
    }

    auto nTrailingBits = dim % BLOCK_SIZE;
    mBlocks[mNumBlocks - 1] = nTrailingBits == 0 ? block
        : block & ((1ull << nTrailingBits) - 1);
}

VecGF2::VecGF2(const std::vector<bool>& symbols)
    : VecGF2(symbols.size(), 0)
{
    block_t mask = 1;
    mBlocks[0] = 0;

    for (size_t pos = 0, blockIdx = 0; pos < mDimension; pos++) {
        if (pos > 0 && pos % BLOCK_SIZE == 0) {
            mask = 1;
            mBlocks[++blockIdx] = 0;
        }

        if (symbols[pos]) {
            mBlocks[blockIdx] |= mask;
        }

        mask <<= 1;
    }
}

VecGF2::VecGF2(const VecGF2& other)
    : VecGF2(other.mDimension)
{
    std::copy_n(other.mBlocks.get(), mNumBlocks, mBlocks.get());
}

VecGF2::VecGF2(VecGF2&& other) noexcept
    : mNumBlocks(other.mNumBlocks)
    , mDimension(other.mDimension)
    , mBlocks(std::move(other.mBlocks))
{}

VecGF2& VecGF2::operator=(const VecGF2& other)
{
    mNumBlocks = other.mNumBlocks;
    mDimension = other.mDimension;
    mBlocks.reset(new block_t[mNumBlocks]);
    std::copy_n(other.mBlocks.get(), mNumBlocks, mBlocks.get());
    return *this;
}

VecGF2& VecGF2::operator=(VecGF2&& other) noexcept
{
    mNumBlocks = other.mNumBlocks;
    mDimension = other.mDimension;
    mBlocks = std::move(other.mBlocks);
    return *this;
}

size_t VecGF2::Dimension() const
{
    return mDimension;
}

size_t VecGF2::HammingWeight() const
{
    size_t w = 0;
    for (size_t i = 0; i < mDimension; i++) {
        w += this->operator[](i);
    }
    return w;
}

std::optional<size_t> VecGF2::End() const
{
    size_t i = mNumBlocks - 1;
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

    return i * BLOCK_SIZE + highestOnePos - 1;
}

bool VecGF2::operator[](size_t pos) const
{
    auto offset = pos / BLOCK_SIZE;
    auto index = pos % BLOCK_SIZE;
    return (mBlocks[offset] & (1ull << index)) != 0;
}

void VecGF2::Set(size_t pos, bool value)
{
    auto offset = pos / BLOCK_SIZE;
    auto index = pos % BLOCK_SIZE;

    if (value) {
        mBlocks[offset] |= 1ull << index;
    }
    else {
        mBlocks[offset] &= ONES ^ (1ull << index);
    }
}

bool VecGF2::operator*(const VecGF2& other) const
{
    auto max_dim = std::max(mDimension, other.mDimension);
    VecGF2 result(max_dim);

    size_t count = 0;

    for (size_t i = 0; i < result.mNumBlocks; i++) {
        auto thisBlock = mBlocks[i];
        auto otherBlock = other.mBlocks[i];
        count += std::bitset<BLOCK_SIZE>(thisBlock & otherBlock).count();
    }

    return count % 2;
}

VecGF2 VecGF2::operator+(const VecGF2& other) const
{
    auto maxDimension = std::max(mDimension, other.mDimension);
    VecGF2 result(maxDimension);

    for (size_t i = 0; i < result.mNumBlocks; i++) {
        result.mBlocks[i] = mBlocks[i] ^ other.mBlocks[i];
    }

    return result;
}

VecGF2& VecGF2::operator+=(const VecGF2& other)
{
    for (size_t i = 0; i < other.mNumBlocks; i++) {
        mBlocks[i] ^= other.mBlocks[i];
    }

    return *this;
}

bool VecGF2::operator==(const VecGF2& other) const
{
    if (mDimension != other.Dimension()) {
        return false;
    }

    for (size_t i = 0; i < mDimension; i++) {
        if (this->operator[](i) != other[i]) {
            return false;
        }
    }
    return true;
}

bool VecGF2::operator!=(const VecGF2& other) const
{
    return !(*this == other);
}
