#pragma once

#include "Vector.h"

namespace math::geom {

    template <size_t Rows, size_t Columns>
    struct Matrix : public std::array<Vector<Columns>, Rows> {
        using value_type = typename Vector<Columns>::value_type;
        using std::array<Vector<Columns>, Rows>::array;
        Matrix(std::initializer_list<std::initializer_list<value_type>> init);

        Matrix& operator*=(value_type scalar);
        Matrix& operator/=(value_type scalar);

        Matrix& operator+=(const Matrix& other);
        Matrix& operator-=(const Matrix& other);

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

        value_type det() const;
    };

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

    inline Matrix2x2 rotation_matrix_2d(double angle_rad);
    inline Matrix3x3 rotation_matrix_x(double angle_rad);
    inline Matrix3x3 rotation_matrix_y(double angle_rad);
    inline Matrix3x3 rotation_matrix_z(double angle_rad);

    inline Matrix4x4 orthographic_projection(double left, double right, double bottom, double top, double near, double far);
    inline Matrix4x4 perspective_projection(double fov_rad, double aspect, double near, double far);

    inline Matrix2x2 scaling_matrix(double sx, double sy);
    inline Matrix3x3 scaling_matrix(double sx, double sy, double sz);
    inline Matrix4x4 scaling_matrix(double sx, double sy, double sz, double sw);

    inline Matrix4x4 scale_around_point(float sx, float sy, float sz, float cx, float cy, float cz);

    inline Matrix3x3 translation_matrix(double x, double y);
    inline Matrix4x4 translation_matrix(double x, double y, double z);

} // namespace math::geom

#include "Matrix.hpp"
