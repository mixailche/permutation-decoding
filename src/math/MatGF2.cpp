#include <stdexcept>
#include <algorithm>
#include <utility>
#include <queue>

#include "math/MathUtils.h"
#include "math/MatGF2.h"

using math::MatGF2;
using math::VecGF2;

MatGF2::MatGF2(size_t numRows, size_t numCols)
    : mNumRows(numRows)
    , mNumCols(numCols)
    , mRows(numRows)
{}

MatGF2 MatGF2::Identity(size_t n)
{
    MatGF2 mat(n, n);
    for (size_t i = 0; i < n; i++) {
        mat[i] = VecGF2(n, 0);
        mat.Set(i, i, 1);
    }
    return mat;
}

MatGF2 MatGF2::Zeros(size_t numRows, size_t numCols)
{
    MatGF2 mat(numRows, numCols);
    std::generate(mat.mRows.begin(), mat.mRows.end(), [=] { return VecGF2(numCols, 0); });
    return mat;
}

MatGF2::MatGF2(const MatGF2& other)
    : mNumRows(other.NumRows())
    , mNumCols(other.NumCols())
    , mRows(other.mRows)
{}

MatGF2::MatGF2(MatGF2 && other) noexcept
    : mNumRows(other.NumRows())
    , mNumCols(other.NumCols())
    , mRows(std::move(other.mRows))
{}

MatGF2::MatGF2(const std::vector<VecGF2>& rows)
    : MatGF2(rows.size(), rows[0].Dimension())
{
    for (size_t i = 0; i < mNumRows; i++) {
        if (rows[i].Dimension() != mNumCols) {
            throw std::invalid_argument("Rows have distinct length");
        }
        mRows[i] = rows[i];
    }
}

MatGF2& MatGF2::operator=(const MatGF2& other)
{
    mNumRows = other.NumRows();
    mNumCols = other.NumCols();
    mRows = other.mRows;
    return *this;
}

MatGF2& MatGF2::operator=(MatGF2&& other) noexcept
{
    mNumRows = other.NumRows();
    mNumCols = other.NumCols();
    mRows = std::move(other.mRows);
    return *this;
}

size_t MatGF2::NumRows() const
{
    return mNumRows;
}

size_t MatGF2::NumCols() const
{
    return mNumCols;
}

const VecGF2& MatGF2::operator[](size_t i) const
{
    return mRows[i];
}

VecGF2& MatGF2::operator[](size_t i)
{
    return mRows[i];
}

bool MatGF2::Determinant() const
{
    if (mNumRows != mNumCols) {
        throw std::invalid_argument("Matrix must be square");
    }

    MatGF2 mat = *this;

    for (size_t i = mNumRows; i-- > 0;) {
        size_t j = i;
        while (!mat[j][i]) {
            if (j-- == 0) {
                return 0;
            }
        }
        if (i != j) {
            std::swap(mat[i], mat[j]);
        }
        for (size_t k = i; k-- > 0;) {
            if (mat[k][i]) {
                mat[k] += mat[i];
            }
        }
    }

    return 1;
}

MatGF2 MatGF2::Invert() const
{
    if (mNumRows != mNumCols) {
        throw std::invalid_argument("Matrix must be square");
    }

    auto dim = mNumRows;
    auto inv = Identity(dim);
    MatGF2 mat = *this;

    for (size_t i = dim; i-- > 0;) {
        // move row vector with 1 at ith position to ith row
        size_t j = i;
        while (!mat[j][i]) {
            if (j-- == 0) {
                throw std::invalid_argument("Singular matrix");
            }
        }
        if (i != j) {
            std::swap(mat[i], mat[j]);
            std::swap(inv[i], inv[j]);
        }

        // make upper row vectors have 0 at ith position
        for (size_t k = i; k-- > 0;) {
            if (mat[k][i]) {
                mat[k] += mat[i];
                inv[k] += inv[i];
            }
        }
    }

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = i + 1; j < dim; j++) {
            if (mat[j][i]) {
                mat[j] += mat[i];
                inv[j] += inv[i];
            }
        }
    }

    return inv;
}

MatGF2 MatGF2::Transpose() const
{
    auto transposed = Zeros(mNumCols, mNumRows);
    for (size_t i = 0; i < transposed.NumRows(); i++) {
        for (size_t j = 0; j < transposed.NumCols(); j++) {
            transposed.Set(i, j, mRows[j][i]);
        }
    }
    return transposed;
}

void MatGF2::SelfTranspose()
{
    for (size_t i = 0; i < mNumRows; i++) {
        for (size_t j = 0; j < mNumCols; j++) {
            auto tmp = mRows[i][j];
            Set(i, j, mRows[j][i]);
            Set(j, i, tmp);
        }
    }
}

VecGF2 math::MatGF2::operator*(const VecGF2& vec) const
{
    if (vec.Dimension() != mNumCols) {
        throw std::invalid_argument("Incorrect matrxi and vector dimensions");
    }

    std::vector<bool> result(mNumCols);
    for (size_t i = 0; i < mNumRows; i++) {
        result[i] = mRows[i] * vec;
    }
    return VecGF2(result);
}

MatGF2 MatGF2::operator*(const MatGF2& other) const
{
    if (other.NumRows() != mNumCols) {
        throw std::invalid_argument("Incorrect multipliers dimensions");
    }

    auto result = Zeros(mNumRows, other.NumCols());
    auto otherT = other.Transpose();

    for (size_t i = 0; i < result.NumRows(); i++) {
        for (size_t j = 0; j < result.NumCols(); j++) {
            result.Set(i, j, mRows[i] * otherT[j]);
        }
    }

    return result;
}

MatGF2& MatGF2::operator*=(const MatGF2& other)
{
    return (*this = *this * other);
}

bool MatGF2::operator==(const MatGF2& other) const
{
    if (mNumRows != other.NumRows() || mNumCols != other.NumCols()) {
        return false;
    }

    for (size_t i = 0; i < mNumRows; i++) {
        if (mRows[i] != other[i]) {
            return false;
        }
    }
    return true;
}

bool MatGF2::operator!=(const MatGF2& other) const
{
    return !(*this == other);
}

void MatGF2::Set(size_t i, size_t j, bool value)
{
    mRows[i].Set(j, value);
}

void MatGF2::AppendRow(VecGF2 row)
{
    mRows.push_back(std::move(row));
}

void MatGF2::RunGaussianElimination()
{
    std::vector<std::queue<size_t>> buckets(mNumCols);

    for (size_t r = 0; r < mNumRows; r++) {
        if (auto c = mRows[r].End(); c.has_value()) {
            buckets[*c].push(r);
        }
    }

    for (size_t c = mNumCols; c-- > 0;) {
        auto& bucket = buckets[c];
        if (bucket.empty()) {
            continue;
        }

        auto r = bucket.front(); bucket.pop();
        auto& base = mRows[r];

        while (!bucket.empty()) {
            auto r1 = bucket.front(); bucket.pop();
            auto& row = mRows[r1];
            row += base;
            if (auto c1 = row.End(); c1.has_value()) {
                buckets[*c1].push(r1);
            }
        }
    }
}


