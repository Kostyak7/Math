#pragma once 

#include "IMatrixFromVector.h"

namespace math::linal {

    /**
     @brief Ленточная матрица особого вида
     @details Симметричная, хранятся только строки начиная с диагонального элемента
    */
    class MATH_EXPORT BandMatrix final: public IMatrixFromVector {
    public:
        BandMatrix() noexcept;
        BandMatrix(size_t size, size_t isl = 1, const value_type& default_value = {});
        BandMatrix(const BandMatrix& matrix);
        BandMatrix(BandMatrix&& matrix) noexcept;
        BandMatrix(std::initializer_list<value_type> list);
        BandMatrix(std::initializer_list<std::initializer_list<value_type>> list);

        BandMatrix& operator=(const BandMatrix& matrix);
        BandMatrix& operator=(BandMatrix&& matrix) noexcept;

        BandMatrix& operator*=(value_type scalar);
        BandMatrix& operator/=(value_type scalar);

        BandMatrix& operator+=(const BandMatrix& other);
        BandMatrix& operator-=(const BandMatrix& other);

        void swap(BandMatrix& matrix) noexcept;

        void reshape(size_t size, size_t isl) override;

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

        size_t get_isl() const;
        void set_isl(size_t width);        
        void shrink();

        value_type det() const override;

        std::vector<std::pair<size_t, complex_value_type>> get_eigenvalues() const override;
        std::vector<std::pair<complex_value_type, std::vector<DVector>>> get_eigenvectors() const override;

        static BandMatrix identity_matrix(size_t n, size_t isl = 1);
        static BandMatrix elementary_matrix_unit(size_t n, size_t i, size_t j, value_type value = 1.0);

    private:
        ConstProxyVector get_row(size_t row) const override;
        ProxyVector get_row(size_t row) override;
        bool check_row_index(size_t row) const override;
        size_t calculate_index(size_t row, size_t col) const override;

    private:
        size_t m_size = 0;
        size_t m_isl = 0;
    };

    void MATH_EXPORT swap(BandMatrix& m1, BandMatrix& m2) noexcept;

    MATH_EXPORT bool operator==(const BandMatrix& m1, const BandMatrix& m2);
    MATH_EXPORT bool operator!=(const BandMatrix& m1, const BandMatrix& m2);

    MATH_EXPORT BandMatrix operator*(const BandMatrix& matrix, BandMatrix::value_type scalar);
    MATH_EXPORT BandMatrix operator*(BandMatrix::value_type scalar, const BandMatrix& matrix);
    MATH_EXPORT BandMatrix operator/(const BandMatrix& matrix, BandMatrix::value_type scalar);

    MATH_EXPORT DVector operator*(const BandMatrix& matrix, const DVector& vector);
    MATH_EXPORT DVector operator*(const DVector& vector, const BandMatrix& matrix);

    MATH_EXPORT BandMatrix operator*(const BandMatrix& m1, const BandMatrix& m2);

    MATH_EXPORT BandMatrix operator+(const BandMatrix& m1, const BandMatrix& m2);
    MATH_EXPORT BandMatrix operator-(const BandMatrix& m1, const BandMatrix& m2);

    MATH_EXPORT BandMatrix operator-(const BandMatrix& matrix);

    BandMatrix MATH_EXPORT inversed(const BandMatrix& matrix);

} // namespace math::linal
