#pragma once 

#include "IMatrixFromVector.h"

namespace math::linal {

    /**
     @brief Обычная плотная матрица 
    */
    class MATH_EXPORT DenseMatrix final: public IMatrixFromVector {
    public:
        DenseMatrix() noexcept;
        DenseMatrix(size_t height, size_t width, const value_type& default_value = {});
        DenseMatrix(const DenseMatrix& matrix);
        DenseMatrix(DenseMatrix&& matrix) noexcept;
        DenseMatrix(std::initializer_list<value_type> list);
        DenseMatrix(std::initializer_list<std::initializer_list<value_type>> list);

        DenseMatrix& operator=(const DenseMatrix& matrix);
        DenseMatrix& operator=(DenseMatrix&& matrix) noexcept;

        DenseMatrix& operator*=(value_type scalar);
        DenseMatrix& operator/=(value_type scalar);

        DenseMatrix& operator+=(const DenseMatrix& other);
        DenseMatrix& operator-=(const DenseMatrix& other);     

        void swap(DenseMatrix& matrix) noexcept;

        void reshape(size_t height, size_t width) override;

        bool is_positive_definite() const override;

        value_type get(size_t row, size_t col) const override;
        void set(size_t row, size_t col, value_type value) override;

        size_t get_width() const override;
        size_t get_height() const override;

        value_type det() const override;

        std::vector<std::pair<size_t, complex_value_type>> get_eigenvalues() const override;
        std::vector<std::pair<complex_value_type, std::vector<DVector>>> get_eigenvectors() const override;

        void swap_rows(size_t r1, size_t r2);
        void swap_columns(size_t c1, size_t c2);

        static DenseMatrix identity_matrix(size_t n);
        static DenseMatrix elementary_matrix_unit(size_t n, size_t m, size_t i, size_t j, double value = 1.0);

    private:
        ConstProxyVector get_row(size_t row) const override;
        ProxyVector get_row(size_t row) override;
        bool check_row_index(size_t row) const override;
        size_t calculate_index(size_t row, size_t col) const override;

    private:
        size_t m_height;
        size_t m_width;
    };

    void MATH_EXPORT swap(DenseMatrix& m1, DenseMatrix& m2) noexcept;

    MATH_EXPORT bool operator==(const DenseMatrix& m1, const DenseMatrix& m2);
    MATH_EXPORT bool operator!=(const DenseMatrix& m1, const DenseMatrix& m2);

    MATH_EXPORT DenseMatrix operator*(const DenseMatrix& matrix, DenseMatrix::value_type scalar);
    MATH_EXPORT DenseMatrix operator*(DenseMatrix::value_type scalar, const DenseMatrix& matrix);
    MATH_EXPORT DenseMatrix operator/(const DenseMatrix& matrix, DenseMatrix::value_type scalar);

    MATH_EXPORT DVector operator*(const DenseMatrix& matrix, const DVector& vector);
    MATH_EXPORT DVector operator*(const DVector& vector, const DenseMatrix& matrix);

    MATH_EXPORT DenseMatrix operator*(const DenseMatrix& m1, const DenseMatrix& m2);

    MATH_EXPORT DenseMatrix operator+(const DenseMatrix& m1, const DenseMatrix& m2);
    MATH_EXPORT DenseMatrix operator-(const DenseMatrix& m1, const DenseMatrix& m2);

    MATH_EXPORT DenseMatrix operator-(const DenseMatrix& matrix);

    DenseMatrix MATH_EXPORT inversed(const DenseMatrix& matrix);

    DenseMatrix MATH_EXPORT transposed(const DenseMatrix& matrix);

} // namespace math::linal
