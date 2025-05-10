#pragma once

#include <vector>
#include "VecGF2.h"

namespace math {
    class MatGF2 {
    public:
        static MatGF2 Identity(size_t n);
        static MatGF2 Zeros(size_t numRows, size_t numCols);

        MatGF2() = default;
        MatGF2(const MatGF2& other);
        MatGF2(MatGF2&& other) noexcept;
        MatGF2(const std::vector<VecGF2>& rows);

        MatGF2& operator=(const MatGF2& other);
        MatGF2& operator=(MatGF2&& other) noexcept;

        ~MatGF2() = default;

        size_t NumRows() const;
        size_t NumCols() const;

        const VecGF2& operator[](size_t i) const;
        VecGF2& operator[](size_t i);

        bool Determinant() const;
        MatGF2 Invert() const;
        MatGF2 Transpose() const;
        void SelfTranspose();

        VecGF2 operator*(const VecGF2& vec) const;
        MatGF2 operator*(const MatGF2& other) const;
        MatGF2& operator*=(const MatGF2& other);

        bool operator==(const MatGF2& other) const;
        bool operator!=(const MatGF2& other) const;

        void Set(size_t i, size_t j, bool value);
        void AppendRow(VecGF2 row);
        void RunGaussianElimination();

    private:
        size_t mNumRows = 0;
        size_t mNumCols = 0;
        std::vector<VecGF2> mRows;

        MatGF2(size_t numRows, size_t numCols);
    };
}