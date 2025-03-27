#pragma once

#include <vector>
#include "Linalg/BinVector.h"

namespace Linalg {
    class BinMatrix {
    public:
        BinMatrix(const std::vector<std::vector<bool>>& symbols);
        ~BinMatrix();

        BinMatrix(BinMatrix&& other) noexcept;
        BinMatrix(const BinMatrix& other);

        BinMatrix operator*(const BinMatrix& other) const;

        void RunGaussianElimination();

        std::vector<bool> SCanColumn(size_t c) const;
        const BinVector& operator[](size_t r) const;

        size_t NRows() const;
        size_t NCols() const;

    private:
        BinMatrix(size_t nRows, size_t nCols);

        size_t mNRows;
        size_t mNCols;
        std::unique_ptr<BinVector[]> mRows;
    };
}
