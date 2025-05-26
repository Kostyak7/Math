template <size_t N>
math::geom::Point<N> math::geom::project_onto(const Point<N>& point, const Line<N>& line) {
    // ...
    return {};
}

template <size_t N>
double math::geom::angle_between(const Vector<N>& v1, const Vector<N>& v2) {
    double norm_product = v1.length() * v2.length();
    if (dcmp(norm_product) == 0) 
        return 0.0;
    return std::acos(clamp(v1.dot(v2) / norm_product, -1.0, 1.0));
}

template <size_t N>
double math::geom::distance_to_line_nd(const Point<N>& point, const Line<N>& line) {
    Vector<N> ab = line.get_vector();
    Vector<N> ap = point - line.p1;
    Vector<N> proj = ab * (ap.dot(ab) / ab.dot(ab));
    Vector<N> perpendicular = ap - proj;
    return perpendicular.length();
}

template <size_t N>
bool math::geom::is_point_on_line(const Point<N>& p, const Line<N>& line, double eps) {
    return distance_to_line(point, line) < eps;
}

template <size_t N>
math::geom::Point<N> math::geom::closest_point_on_line(const Point<N>& point, const Line<N>& line) {
    // ...
    return {};
}
