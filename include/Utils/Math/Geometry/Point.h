#pragma once

#include <Utils/Math/Common.h>

#include <array>

namespace math::geom {

    template <size_t N, typename T = double>
    class Point {
    public:
        using value_type = T;
        using base_container_type = std::array<T, N>;
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
        Point(value_type default_value = 0.0);
        Point(std::initializer_list<value_type> init);
        Point(const base_container_type& container);
        Point(const Point& other);
        template <class T2> Point(const Point<N, T2>& other);
        Point(Point&& other) noexcept;
        ~Point() = default;

        Point& operator=(const Point& other);
        template <class T2> Point& operator=(const Point<N, T2>& other);
        Point& operator=(Point&& other) noexcept;

        Point& operator*=(value_type scalar);
        Point& operator/=(value_type scalar);

        Point& operator+=(const Point& other);
        Point& operator-=(const Point& other);

        void swap(Point& other) noexcept;

        value_type distance_to(const Point& other);

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

    private:
        base_container_type m_data;
    };

    template <size_t N, typename T>
    void swap(Point<N, T>& p1, Point<N, T>& p2) noexcept;

    template <size_t N, typename T>
    Point<N, T> operator*(Point<N, T> point, T scalar);
    template <size_t N, typename T>
    Point<N, T> operator*(T scalar, Point<N, T> point);
    template <size_t N, typename T>
    Point<N, T> operator/(Point<N, T> point, T scalar);

    template <size_t N, typename T>
    Point<N, T> operator+(Point<N, T> p1, const Point<N, T>& p2);
    template <size_t N, typename T>
    Point<N, T> operator-(Point<N, T> p1, const Point<N, T>& p2);

    template <size_t N, typename T>
    Point<N, T> operator-(const Point<N, T>& point);

    template <size_t N, typename T>
    Point<N, T> zero_point();
    template <size_t N, typename T>
    Point<N, T> unit_point(size_t i, T value = 1.0);

    using Point2D = Point<2, double>;
    using Point3D = Point<3, double>;
    using Point4D = Point<4, double>;

} // namespace math::geom

#include "Point.tpp"
