#pragma once 

#include "IMatrix.h"

#include <unordered_map>
#include <utility>
#include <functional>

namespace std {
    template<>
    struct hash<std::pair<size_t, size_t>> {
        size_t operator()(const std::pair<size_t, size_t>& p) const noexcept {
            return std::hash<size_t>{}(p.first) ^ (std::hash<size_t>{}(p.second) << 1);
        }
    };
}

namespace math::linal {

    class SparseMatrix final: public IMatrix {
    public:
        SparseMatrix() noexcept = default;
        SparseMatrix(size_t width, size_t height);
        SparseMatrix(const SparseMatrix& matrix);
        SparseMatrix(SparseMatrix&& matrix) noexcept;

        SparseMatrix& operator=(const SparseMatrix& other);
        SparseMatrix& operator=(SparseMatrix&& other);

        SparseMatrix& operator*=(value_type scalar) noexcept;
        SparseMatrix& operator/=(value_type scalar);

        SparseMatrix& operator+=(const SparseMatrix& other);
        SparseMatrix& operator-=(const SparseMatrix& other);

        bool is_zero() const override;
        bool is_identity() const override;
        bool is_diagonal() const override;
        bool is_symmetrical() const override;
        bool is_upper_triangular() const override;
        bool is_lower_triangular() const override;
        bool is_trapezoidal() const override;

        value_type get(size_t row, size_t col) const override;
        void set(size_t row, size_t col, value_type value) override;

        size_t get_width() const override;
        size_t get_height() const override;

        void shrink();
        void clear() override;

        value_type det() const override;

        std::vector<std::pair<size_t, complex_value_type>> get_eigenvalues() const override;
        std::vector<std::pair<complex_value_type, std::vector<FVector>>> get_eigenvectors() const override;

        void swap(SparseMatrix& other) noexcept;

        static SparseMatrix identity_matrix(size_t n);

        const std::unordered_map<std::pair<size_t, size_t>, value_type>& data() const;
        std::unordered_map<std::pair<size_t, size_t>, value_type> data();

    private:
        size_t m_width = 0;
        size_t m_height = 0;
        std::unordered_map<std::pair<size_t, size_t>, value_type> m_data;

    private:
        friend bool operator==(const SparseMatrix& m1, const SparseMatrix& m2);

        friend FVector operator*(const SparseMatrix& matrix, const FVector& vector);
        friend FVector operator*(const FVector& vector, const SparseMatrix& matrix);

        friend SparseMatrix operator*(const SparseMatrix& m1, const SparseMatrix& m2);

        friend SparseMatrix transposed(const SparseMatrix& matrix);
    };

    void swap(SparseMatrix& m1, SparseMatrix& m2) noexcept;

    bool operator==(const SparseMatrix& m1, const SparseMatrix& m2);
    bool operator!=(const SparseMatrix& m1, const SparseMatrix& m2);

    SparseMatrix operator*(const SparseMatrix& matrix, SparseMatrix::value_type scalar);
    SparseMatrix operator*(SparseMatrix::value_type scalar, const SparseMatrix& matrix);
    SparseMatrix operator/(const SparseMatrix& matrix, SparseMatrix::value_type scalar);

    FVector operator*(const SparseMatrix& matrix, const FVector& vector);
    FVector operator*(const FVector& vector, const SparseMatrix& matrix);

    SparseMatrix operator*(const SparseMatrix& m1, const SparseMatrix& m2);

    SparseMatrix operator+(const SparseMatrix& m1, const SparseMatrix& m2);
    SparseMatrix operator-(const SparseMatrix& m1, const SparseMatrix& m2);

    SparseMatrix operator-(const SparseMatrix& matrix);

    SparseMatrix inversed(const SparseMatrix& matrix);

    SparseMatrix transposed(const SparseMatrix& matrix);

} // namespace math::linal