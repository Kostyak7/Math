template <size_t N>
math::geom::Line<N>::Line(const Vector<N>& vec) {
    p1.fill(0.0);
    std::copy(vec.begin(), vec.end(), p2.begin());
}

template <size_t N>
math::geom::Line<N>::Line(const Point<N>& point) 
    : p2(p2)
{
    p1.fill(0.0);
}

template <size_t N>
math::geom::Line<N>::Line(const Point<N>& p1, const Point<N>& p2) 
    : p1(p1)
    , p2(p2)
{
}

template <size_t N>
typename math::geom::Line<N>::value_type math::geom::Line<N>::length() const {
    return p1.distance_to(p2);
}

template <size_t N>
bool math::geom::Line<N>::contains(const Point<N>& point) const {
    return parallel(Line<N>(p1, point));
}

template <size_t N>
bool math::geom::Line<N>::parallel(const Line<N>& other) const {
    if constexpr (N == 2) {
        Vector2D v1(p1, p2);
        Vector2D v2(other.p1, other.p2);
        return fabs(v1.pseudodot(v2)) < std::numeric_limits<Point2D::value_type>::epsilon();
    }
    else if constexpr (N == 3) {
        Vector3D v1(p1, p2);
        Vector3D v2(other.p1, other.p2);
        Vector3D cross = v1.cross(v2);
        return fabs(cross.x) < std::numeric_limits<Point3D::value_type>::epsilon() &&
            fabs(cross.y) < std::numeric_limits<Point3D::value_type>::epsilon() &&
            fabs(cross.z) < std::numeric_limits<Point3D::value_type>::epsilon();
    }
    return false;
}

template <size_t N>
bool math::geom::Line<N>::collinear(const Line<N>& other) const {
    return parallel(other) && contains(other.p1);
}

template <size_t N>
bool math::geom::Line<N>::perpendicular(const Line<N>& other) const {
    Vector<N> v1(p1, p2);
    Vector<N> v2(other.p1, other.p2);

    return fabs(v1.dot(v2)) < std::numeric_limits<Point3D::value_type>::epsilon();
}

template <size_t N>
math::geom::Vector<N> math::geom::Line<N>::get_vector() const {
    return { p2 - p1 };
}
