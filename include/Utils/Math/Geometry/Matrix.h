#pragma once

#include "Vector.h"

namespace math::geom {

    template <size_t Rows, size_t Columns>
    class Matrix {
    public:
        using value_type = double;
        
        class ConstProxyVector;
        class ProxyVector;

    public:
        Matrix();
        Matrix(const Matrix& matrix);
        Matrix(Matrix&& matrix) noexcept;
        Matrix(std::initializer_list<std::initializer_list<value_type>> init);

        Matrix& operator=(const Matrix& matrix);
        Matrix& operator=(Matrix&& matrix) noexcept;

        Matrix& operator*=(value_type scalar);
        Matrix& operator/=(value_type scalar);

        Matrix& operator+=(const Matrix& other);
        Matrix& operator-=(const Matrix& other);

        void swap(const Matrix& other) noexcept;

        bool is_empty() const;

        bool is_square() const;
        bool is_zero() const;
        bool is_identity() const;
        bool is_diagonal() const;
        bool is_symmetrical() const;
        bool is_upper_triangular() const;
        bool is_lower_triangular() const;
        bool is_trapezoidal() const;
        bool is_positive_definite() const;
        bool is_negative_definite() const;

        size_t get_width() const;
        size_t get_height() const;

        value_type get(size_t row, size_t col) const;
        void set(size_t row, size_t col, value_type value);

        ConstProxyVector operator[](size_t row) const;
        ProxyVector operator[](size_t row);

        void swap_rows(size_t r1, size_t r2);
        void swap_columns(size_t c1, size_t c2);

        value_type det() const;

    private:
        const value_type& _(size_t row, size_t col) const; //direct access
        value_type& _(size_t row, size_t col); // Its not an operator() because in methods it is perceived as an operator,

    private:
        std::array<value_type, Rows * Columns> m_data;
    };

    template <size_t Rows, size_t Columns>
    void swap(Matrix<Rows, Columns>& m1, Matrix<Rows, Columns>& m2) noexcept;

    template <size_t Rows, size_t Columns>
    bool operator==(const Matrix<Rows, Columns>& m1, const Matrix<Rows, Columns>& m2);
    template <size_t Rows, size_t Columns>
    bool operator!=(const Matrix<Rows, Columns>& m1, const Matrix<Rows, Columns>& m2);

    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator*(const Matrix<Rows, Columns>& matrix, typename Matrix<Rows, Columns>::value_type scalar);
    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator*(typename Matrix<Rows, Columns>::value_type scalar, const Matrix<Rows, Columns>& matrix);
    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator/(const Matrix<Rows, Columns>& matrix, typename Matrix<Rows, Columns>::value_type scalar);

    template <size_t Rows, size_t Columns>
    Vector<Rows> operator*(const Matrix<Rows, Columns>& matrix, const Vector<Columns>& vector);
    template <size_t Rows, size_t Columns>
    Vector<Columns> operator*(const Vector<Rows>& vector, const Matrix<Rows, Columns>& matrix);

    template <size_t N, size_t M, size_t K>
    Matrix<N, K> operator*(const Matrix<N, M>& m1, const Matrix<M, K>& m2);

    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator+(const Matrix<Rows, Columns>& m1, const Matrix<Rows, Columns>& m2);
    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator-(const Matrix<Rows, Columns>& m1, const Matrix<Rows, Columns>& m2);

    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator-(const Matrix<Rows, Columns>& matrix);

    template <size_t N>
    Matrix<N, N> inversed(const Matrix<N, N>& matrix);

    template <size_t Rows, size_t Columns>
    Matrix<Columns, Rows> transposed(const Matrix<Rows, Columns>& matrix);

    template <size_t N>
    Matrix<N, N> identity_matrix();
    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> elementary_matrix_unit(size_t i, size_t j, double value = 1.0);

    using Matrix2x2 = Matrix<2, 2>;
    using Matrix3x3 = Matrix<3, 3>;
    using Matrix4x4 = Matrix<4, 4>;

    template <size_t Rows, size_t Columns>
    class Matrix<Rows, Columns>::ConstProxyVector {
    public:
        using const_iterator = typename std::array<value_type, Columns>::const_iterator;

        const_iterator begin() const;
        const_iterator end() const;

    public:
        ConstProxyVector(const_iterator front);

        value_type operator[](size_t col) const;

        size_t size() const;
        bool empty() const;

        operator Vector<Columns>() const;

    private:
        const_iterator m_front;
    };

    template <size_t Rows, size_t Columns>
    bool operator==(const Matrix<Rows, Columns>::ConstProxyVector& v1, const Matrix<Rows, Columns>::ConstProxyVector& v2);
    template <size_t Rows, size_t Columns>
    bool operator!=(const Matrix<Rows, Columns>::ConstProxyVector& v1, const Matrix<Rows, Columns>::ConstProxyVector& v2);

    template <size_t Rows, size_t Columns>
    class Matrix<Rows, Columns>::ProxyVector {
    public:
        using const_iterator = typename std::array<value_type, Columns>::const_iterator;
        using iterator = typename std::vector<value_type, Columns>::iterator;

        const_iterator begin() const;
        const_iterator end() const;
        iterator begin();
        iterator end();

    public:
        ProxyVector(iterator front);

        value_type operator[](size_t col) const;
        value_type& operator[](size_t col);

        size_t size() const;
        bool empty() const;

        operator Vector<Columns>() const;

    private:
        iterator m_front;
    };

    template <size_t Rows, size_t Columns>
    bool operator==(const Matrix<Rows, Columns>::ProxyVector& v1, const Matrix<Rows, Columns>::ProxyVector& v2);
    template <size_t Rows, size_t Columns>
    bool operator!=(const Matrix<Rows, Columns>::ProxyVector& v1, const Matrix<Rows, Columns>::ProxyVector& v2);

} // namespace math::geom

#include "Matrix.tpp"
