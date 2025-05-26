#pragma once

#include <Vector.h>

namespace math::geom {

    template <size_t N>
    class Polygon {
    public:
        using value_type = typename Point<N>::value_type;

        Polygon() = default;
        Polygon(const std::vector<Point<N>>& points)
            : m_points(points) 
        {
        }

        const std::vector<Point<N>>& vertices() const {
            return m_points;
        }

        std::vector<Point<N>>& vertices() {
            return m_points;
        }

        size_t size() const { 
            return m_points.size(); 
        }

        bool is_planar(double eps = 1e-9) const {
            if (m_points.size() < 4)
                return true;
            // ...
            return true;
        }

        value_type area() const {
            // ...
            return {};
        }

    private:
        std::vector<Point<N>> m_points;
    };

    using Polygon2D = Polygon<2>;
    using Polygon3D = Polygon<3>;

    inline Vector3D normal_vector(const Polygon3D& polygon) const {
        if (polygon.size() < 3)
            return zero_vector<3>();
        Vector3D normal = zero_vector<3>();
        auto points = polygon.vertices();
        for (size_t i = 0; i < points.size(); ++i) {
            const Point3D& current = points[i];
            const Point3D& next = points[(i + 1) % polygon.size()];
            normal += cross_product(current, next);
        }
        return normalized(normal);
    }

} // namespace math::geom

//#include "Polygon.tpp"
