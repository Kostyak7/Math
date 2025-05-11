#pragma once

#include "Vector.h"

namespace math::geom {

    template <size_t N>
    struct Line {
        using value_type = typename Point<N>::value_type;
        Point<N> p1, p2;

        Line(const Vector<N>& vec);
        Line(const Point<N>& point);
        Line(const Point<N>& p1, const Point<N>& p2);

        value_type length() const;

        bool contains(const Point<N>& point) const;
        bool parallel(const Line<N>& other) const;
        bool collinear(const Line<N>& other) const;
        bool perpendicular(const Line<N>& other) const;

        Vector<N> get_vector() const;
    };

    using Line2D = Line<2>;
    using Line3D = Line<3>;

} // namespace math::geom

#include "Line.hpp"
