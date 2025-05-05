#pragma once 

#include "Vector.h"

namespace math::geom {

    class Plane {
    public:
        using value_type = Point3D::value_type;

        value_type a, b, c, d; // ax + by + cz + d = 0

        Plane(const Point3D& p1, const Point3D& p2, const Point3D& p3);

        Vector3D normal_vector() const;
        value_type distance_to(const Point3D& point) const;        
        bool contains(const Point3D& point) const;
        bool contains(const Vector3D& vec) const;
    };

}