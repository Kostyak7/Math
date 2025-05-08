#pragma once 

#include "IMatrix.h"

namespace math::linal {

    /**
     @brief Обычная плотная матрица 
    */
    class DenseMatrix final: public IMatrix,
                             public std::vector<FVector> {
    public:
        using value_type = IMatrix::value_type;
        using std::vector<FVector>::vector;
        DenseMatrix(size_t height, size_t width, const value_type& default_value = {});

    public:
        DenseMatrix& operator*=(value_type scalar) noexcept;
        DenseMatrix& operator/=(value_type scalar);

        DenseMatrix& operator+=(const DenseMatrix& other);
        DenseMatrix& operator-=(const DenseMatrix& other);     

        bool is_positive_definite() const override;

        value_type get(size_t row, size_t col) const override;
        void set(size_t row, size_t col, value_type value) override;

        size_t get_width() const override;
        size_t get_height() const override;

        value_type det() const override;

        std::vector<std::pair<size_t, complex_value_type>> get_eigenvalues() const override;
        std::vector<std::pair<complex_value_type, std::vector<FVector>>> get_eigenvectors() const override;

        static DenseMatrix identity_matrix(size_t n);
        static DenseMatrix elementary_matrix_unit(size_t n, size_t m, size_t i, size_t j);
    };

    bool operator==(const DenseMatrix& m1, const DenseMatrix& m2);
    bool operator!=(const DenseMatrix& m1, const DenseMatrix& m2);

    DenseMatrix operator*(const DenseMatrix& matrix, DenseMatrix::value_type scalar);
    DenseMatrix operator*(DenseMatrix::value_type scalar, const DenseMatrix& matrix);
    DenseMatrix operator/(const DenseMatrix& matrix, DenseMatrix::value_type scalar);

    FVector operator*(const DenseMatrix& matrix, const FVector& vector);
    FVector operator*(const FVector& vector, const DenseMatrix& matrix);

    DenseMatrix operator*(const DenseMatrix& m1, const DenseMatrix& m2);

    DenseMatrix operator+(const DenseMatrix& m1, const DenseMatrix& m2);
    DenseMatrix operator-(const DenseMatrix& m1, const DenseMatrix& m2);

    DenseMatrix operator-(const DenseMatrix& matrix);

    DenseMatrix inversed(const DenseMatrix& matrix);

    DenseMatrix transposed(const DenseMatrix& matrix);

} // namespace math::linal