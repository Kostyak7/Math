#pragma once

#include "Vector.h"

namespace math::geom {

    template <size_t Rows, size_t Columns>
    class Matrix {
    public:
        using value_type = double;
        using base_container_type = std::array<value_type, Rows* Columns>;

        class ConstRowIterator;
        class RowIterator;

        using const_iterator = ConstRowIterator;
        using iterator = RowIterator;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;
        using reverse_iterator = std::reverse_iterator<iterator>;

        iterator begin() noexcept;
        const_iterator begin() const noexcept;
        iterator end() noexcept;
        const_iterator end() const noexcept;
        reverse_iterator rbegin() noexcept;
        const_reverse_iterator rbegin() const noexcept;
        reverse_iterator rend() noexcept;
        const_reverse_iterator rend() const noexcept;
        const_iterator cbegin() const noexcept;
        const_iterator cend() const noexcept;
        const_reverse_iterator crbegin() const noexcept;
        const_reverse_iterator crend() const noexcept;

        class ConstProxyVector;
        class ProxyVector;

        ConstProxyVector front() const noexcept;
        ProxyVector front() noexcept;
        ConstProxyVector back() const noexcept;
        ProxyVector back() noexcept;

    public:
        Matrix() = default;
        Matrix(std::initializer_list<std::initializer_list<value_type>> init);
        Matrix(const base_container_type& container);
        Matrix(const Matrix& other);
        Matrix(Matrix&& other) noexcept;
        ~Matrix() = default;

        Matrix& operator=(const Matrix& other);
        Matrix& operator=(Matrix&& other) noexcept;

        Matrix& operator*=(value_type scalar);
        Matrix& operator/=(value_type scalar);

        Matrix& operator+=(const Matrix& other);
        Matrix& operator-=(const Matrix& other);

        void swap(Matrix& other) noexcept;

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
        base_container_type m_data;
    };

    template <size_t Rows, size_t Columns>
    void swap(Matrix<Rows, Columns>& m1, Matrix<Rows, Columns>& m2) noexcept;

    template <size_t Rows, size_t Columns>
    bool operator==(const Matrix<Rows, Columns>& m1, const Matrix<Rows, Columns>& m2);
    template <size_t Rows, size_t Columns>
    bool operator!=(const Matrix<Rows, Columns>& m1, const Matrix<Rows, Columns>& m2);

    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator*(Matrix<Rows, Columns> matrix, typename Matrix<Rows, Columns>::value_type scalar);
    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator*(typename Matrix<Rows, Columns>::value_type scalar, Matrix<Rows, Columns> matrix);
    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator/(Matrix<Rows, Columns> matrix, typename Matrix<Rows, Columns>::value_type scalar);

    template <size_t Rows, size_t Columns>
    Vector<Rows> operator*(const Matrix<Rows, Columns>& matrix, const Vector<Columns>& vector);
    template <size_t Rows, size_t Columns>
    Vector<Columns> operator*(const Vector<Rows>& vector, const Matrix<Rows, Columns>& matrix);

    template <size_t N, size_t M, size_t K>
    Matrix<N, K> operator*(const Matrix<N, M>& m1, const Matrix<M, K>& m2);

    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator+(Matrix<Rows, Columns> m1, const Matrix<Rows, Columns>& m2);
    template <size_t Rows, size_t Columns>
    Matrix<Rows, Columns> operator-(Matrix<Rows, Columns> m1, const Matrix<Rows, Columns>& m2);

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
    class Matrix<Rows, Columns>::ConstRowIterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = ConstProxyVector;
        using difference_type = std::ptrdiff_t;
        using pointer = const ConstProxyVector*;
        using reference = const ConstProxyVector&;

        ConstRowIterator(const Matrix* matrix, size_t row);

        reference operator*() const;
        pointer operator->() const;

        ConstRowIterator& operator++();
        ConstRowIterator operator++(int);
        ConstRowIterator& operator--();
        ConstRowIterator operator--(int);

        ConstRowIterator operator+(difference_type n) const;
        ConstRowIterator operator-(difference_type n) const;
        template <class Iter> difference_type operator-(const Iter& other) const;

        template <class Iter> bool operator==(const Iter& other) const;
        template <class Iter> bool operator!=(const Iter& other) const;
        template <class Iter> bool operator<(const Iter& other) const;
        template <class Iter> bool operator>(const Iter& other) const;
        template <class Iter> bool operator<=(const Iter& other) const;
        template <class Iter> bool operator>=(const Iter& other) const;

        ConstRowIterator& operator+=(difference_type n);
        ConstRowIterator& operator-=(difference_type n);

        reference operator[](difference_type n) const;

    private:
        const Matrix* m_matrix;
        size_t m_row;
        mutable ConstProxyVector m_proxy;
    };

    template <size_t Rows, size_t Columns>
    class Matrix<Rows, Columns>::RowIterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = ProxyVector;
        using difference_type = std::ptrdiff_t;
        using pointer = ProxyVector*;
        using reference = ProxyVector&;

        RowIterator(Matrix* matrix, size_t row);

        const reference operator*() const;
        reference operator*();

        const pointer operator->() const;
        pointer operator->();

        RowIterator& operator++();
        RowIterator operator++(int);
        RowIterator& operator--();
        RowIterator operator--(int);

        RowIterator operator+(difference_type n) const;
        RowIterator operator-(difference_type n) const;
        template <class Iter> difference_type operator-(const Iter& other) const;

        template <class Iter> bool operator==(const Iter& other) const;
        template <class Iter> bool operator!=(const Iter& other) const;
        template <class Iter> bool operator<(const Iter& other) const;
        template <class Iter> bool operator>(const Iter& other) const;
        template <class Iter> bool operator<=(const Iter& other) const;
        template <class Iter> bool operator>=(const Iter& other) const;

        RowIterator& operator+=(difference_type n);
        RowIterator& operator-=(difference_type n);

        const reference operator[](difference_type n) const;
        reference operator[](difference_type n);

    private:
        Matrix* m_matrix;
        size_t m_row;
        mutable ProxyVector m_proxy;
    };

    template <size_t Rows, size_t Columns>
    class Matrix<Rows, Columns>::ConstProxyVector {
    public:
        using const_iterator = typename base_container_type::const_iterator;
        using const_reverse_iterator = typename base_container_type::const_reverse_iterator;

        const_iterator begin() const noexcept;
        const_iterator end() const noexcept;
        const_reverse_iterator rbegin() const noexcept;
        const_reverse_iterator rend() const noexcept;
        const_iterator cbegin() const noexcept;
        const_iterator cend() const noexcept;
        const_reverse_iterator crbegin() const noexcept;
        const_reverse_iterator crend() const noexcept;

        const value_type& front() const noexcept;
        const value_type& back() const noexcept;

    public:
        ConstProxyVector(const_iterator front);

        size_t size() const;
        bool empty() const;

        bool is_zero() const;

        const value_type& operator[](size_t pos) const noexcept;

        value_type dot(const Vector<Columns>& other) const;
        value_type dot(const ProxyVector& other) const;
        value_type norm() const;
        value_type length() const;

        operator Vector<Columns>() const;

    private:
        const_iterator m_front;
    };

    template <size_t Rows, size_t Columns>
    bool operator==(typename const Matrix<Rows, Columns>::ConstProxyVector& v1, typename const Matrix<Rows, Columns>::ConstProxyVector& v2);
    template <size_t Rows, size_t Columns>
    bool operator!=(typename const Matrix<Rows, Columns>::ConstProxyVector& v1, typename const Matrix<Rows, Columns>::ConstProxyVector& v2);

    template <size_t Rows, size_t Columns>
    class Matrix<Rows, Columns>::ProxyVector {
    public:
        using const_iterator = typename base_container_type::const_iterator;
        using iterator = typename base_container_type::iterator;
        using const_reverse_iterator = typename base_container_type::const_reverse_iterator;
        using reverse_iterator = typename base_container_type::reverse_iterator;

        iterator begin() noexcept;
        const_iterator begin() const noexcept;
        iterator end() noexcept;
        const_iterator end() const noexcept;
        reverse_iterator rbegin() noexcept;
        const_reverse_iterator rbegin() const noexcept;
        reverse_iterator rend() noexcept;
        const_reverse_iterator rend() const noexcept;
        const_iterator cbegin() const noexcept;
        const_iterator cend() const noexcept;
        const_reverse_iterator crbegin() const noexcept;
        const_reverse_iterator crend() const noexcept;

        const value_type& front() const noexcept;
        value_type& front() noexcept;
        const value_type& back() const noexcept;
        value_type& back() noexcept;

    public:
        ProxyVector(iterator front);
        ProxyVector(const ConstProxyVector& proxy);

        ProxyVector& operator*=(value_type scalar);
        ProxyVector& operator/=(value_type scalar);

        ProxyVector& operator+=(const ProxyVector& other);
        ProxyVector& operator-=(const ProxyVector& other);

        size_t size() const;
        bool empty() const;

        void fill(value_type value);

        bool is_zero() const;

        const value_type& operator[](size_t pos) const noexcept;
        value_type& operator[](size_t pos) noexcept;

        value_type dot(const Vector<Columns>& other) const;
        value_type dot(const ProxyVector& other) const;
        value_type norm() const;
        value_type length() const;

        ProxyVector& normalize();

        operator Vector<Columns>() const;

    private:
        iterator m_front;
    };

    template <size_t Rows, size_t Columns>
    bool operator==(typename const Matrix<Rows, Columns>::ProxyVector& v1, typename const Matrix<Rows, Columns>::ProxyVector& v2);
    template <size_t Rows, size_t Columns>
    bool operator!=(typename const Matrix<Rows, Columns>::ProxyVector& v1, typename const Matrix<Rows, Columns>::ProxyVector& v2);

} // namespace math::geom

#include "Matrix.tpp"
