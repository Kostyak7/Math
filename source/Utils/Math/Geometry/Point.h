#pragma once

#include <Utils/Math/Common.h>

#include <array>

namespace math::geom {

    template <size_t N, typename T = double>
    struct Point : public std::array<T, N> {
        using value_type = T;
        using std::array<T, N>::array;

        value_type distance_to(const Point& other);

        value_type x() const;
        value_type y() const;
        value_type z() const;
        value_type w() const;
    };

    template <size_t N, typename T>
    Point<N, T> operator-(const Point<N, T>& point);
    template <size_t N, typename T>
    Point<N, T> operator+(const Point<N, T>& p1, const Point<N, T>& p2);
    template <size_t N, typename T>
    Point<N, T> operator-(const Point<N, T>& p1, const Point<N, T>& p2);
    template <size_t N, typename T>
    Point<N, T> operator*(const Point<N, T>& point, T scalar);
    template <size_t N, typename T>
    Point<N, T> operator*(T scalar, const Point<N, T>& point);
    template <size_t N, typename T>
    Point<N, T> operator/(const Point<N, T>& point, T scalar);


    using Point2D = Point<2, double>;
    using Point3D = Point<3, double>;
    using Point4D = Point<4, double>;

} // namespace math::geom

#include "Point.hpp"
