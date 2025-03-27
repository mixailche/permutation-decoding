#pragma once

#include <vector>
#include <optional>
#include <memory>

namespace Linalg {
    class BinVector {
    public:
        BinVector(size_t dim, bool value);
        BinVector(const std::vector<bool>& symbols);

        BinVector(const BinVector& other);
        BinVector(BinVector&& other) noexcept;

        BinVector& operator=(const BinVector& other);
        BinVector& operator=(BinVector&& other) noexcept;

        size_t Dimension() const;
        std::optional<size_t> end() const;

        bool operator[](size_t pos) const;
        void Set(size_t pos, bool value);

        bool operator*(const BinVector& other) const;

        BinVector operator+(const BinVector& other) const;
        BinVector& operator+=(const BinVector& other);

    private:
        BinVector(size_t dim);

        size_t mDimension;
        size_t mNBlocks;

        using block_t = uint64_t;
        std::unique_ptr<block_t[]> mBlocks;

        constexpr static size_t BlockSize = 8 * sizeof(block_t);
        constexpr static block_t Ones = std::numeric_limits<block_t>::max();
    };
}
