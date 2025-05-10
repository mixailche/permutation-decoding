#pragma once

#include <vector>
#include <optional>
#include <memory>

namespace math {
    class VecGF2 {
    public:
        VecGF2();
        VecGF2(size_t dim, bool value);
        VecGF2(const std::vector<bool>& symbols);

        VecGF2(const VecGF2& other);
        VecGF2(VecGF2&& other) noexcept;

        VecGF2& operator=(const VecGF2& other);
        VecGF2& operator=(VecGF2&& other) noexcept;

        size_t HammingWeight() const;
        size_t Dimension() const;
        std::optional<size_t> End() const;

        bool operator[](size_t pos) const;
        void Set(size_t pos, bool value);

        bool operator*(const VecGF2& other) const;
        
        VecGF2 operator+(const VecGF2& other) const;
        VecGF2& operator+=(const VecGF2& other);

        bool operator==(const VecGF2& other) const;
        bool operator!=(const VecGF2& other) const;

    private:
        size_t mDimension;
        size_t mNumBlocks;

        using block_t = uint64_t;
        std::unique_ptr<block_t[]> mBlocks;

        constexpr static size_t BLOCK_SIZE = 8 /* bits per byte */ * sizeof(block_t);
        constexpr static block_t ONES = std::numeric_limits<block_t>::max();
    
        VecGF2(size_t dimension);
    };
}
