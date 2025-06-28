#pragma once

#include "Point.h"

namespace math::geom {

    template <size_t N>
    class Vector {
    public:
        using value_type = double;
        using base_container_type = std::array<value_type, N>;
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
        Vector(value_type default_value = 0.0);
        Vector(const Point<N, value_type>& to);
        Vector(const Point<N, value_type>& from, const Point<N, value_type>& to);
        Vector(std::initializer_list<value_type> init);
        Vector(const base_container_type& container);
        Vector(const Vector& other);
        Vector(Vector&& other) noexcept;
        ~Vector() = default;

        Vector& operator=(const Vector& other);
        Vector& operator=(Vector&& other) noexcept;

        Vector& operator*=(value_type scalar);
        Vector& operator/=(value_type scalar);

        Vector& operator+=(const Vector& other);
        Vector& operator-=(const Vector& other);

        void swap(Vector& other) noexcept;

        Point<N> get_end_point() const;

        size_t size() const;
        bool empty() const;

        void fill(value_type value);

        bool is_zero() const;

        const value_type& operator[](size_t pos) const noexcept;
        value_type& operator[](size_t pos) noexcept;

        const value_type& at(size_t pos) const;
        value_type& at(size_t pos);

        value_type x() const;
        value_type y() const;
        value_type z() const;
        value_type w() const;

        value_type dot(const Vector& other) const;
        value_type norm() const;
        value_type length() const;

        Vector& normalize();

    private:
        base_container_type m_data;
    };

    template <size_t N>
    void swap(Vector<N>& v1, Vector<N>& v2) noexcept;

    template <size_t N>
    Vector<N> operator*(Vector<N> vec, typename Vector<N>::value_type scalar);
    template <size_t N>
    Vector<N> operator*(typename Vector<N>::value_type scalar, Vector<N> vec);
    template <size_t N>
    Vector<N> operator/(Vector<N> vec, typename Vector<N>::value_type scalar);

    template <size_t N>
    Vector<N> operator+(Vector<N> v1, const Vector<N>& v2);
    template <size_t N>
    Vector<N> operator-(Vector<N> v1, const Vector<N>& v2);

    template <size_t N>
    Vector<N> operator-(const Vector<N>& vec);

    template <size_t N>
    Vector<N> normalized(Vector<N> vector);

    template <size_t N>
    Vector<N> zero_vector();
    template <size_t N>
    Vector<N> unit_vector(size_t i, typename Vector<N>::value_type value = 1.0);

    using Vector2D = Vector<2>;
    using Vector3D = Vector<3>;
    using Vector4D = Vector<4>;

    inline Vector2D::value_type pseudodot(const Vector2D& v1, const Vector2D& v2);

    inline Vector3D cross_product(const Vector3D& v1, const Vector3D& v2);

} // namespace math::geom

#include "Vector.tpp"
