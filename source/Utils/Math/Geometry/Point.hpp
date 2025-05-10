template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type math::geom::Point<N, T>::distance_to(const Point& other) {
    value_type res2 = 0;
    for (size_t i = 0; i < N; ++i) {
        const auto di = (*this)[i] - other[i];
        res2 += di * di;
    }
    return std::sqrt(res2);
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type math::geom::Point<N, T>::x() const {
    if constexpr (N >= 1) {
        return (*this)[0];
    }
    return {};
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type math::geom::Point<N, T>::y() const {
    if constexpr (N >= 2) {
        return (*this)[1];
    }
    return {};
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type math::geom::Point<N, T>::z() const {
    if constexpr (N >= 3) {
        return (*this)[2];
    }
    return {};
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type math::geom::Point<N, T>::w() const {
    if constexpr (N >= 4) {
        return (*this)[3];
    }
    return {};
}

template <size_t N, typename T>
typename math::geom::Point<N, T> math::geom::operator-(const Point<N, T>& point) {
    Point<N, T> res;
    for (size_t i = 0; i < N; ++i) {
        res[i] = - point[i];
    }
    return res;
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::operator+(const Point<N, T>& p1, const Point<N, T>& p2) {
    Point<N, T> res = p1;
    for (size_t i = 0; i < N; ++i) {
        res[i] += p2[i];
    }
    return res;
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::operator-(const Point<N, T>& p1, const Point<N, T>& p2) {
    Point<N, T> res = p1;
    for (size_t i = 0; i < N; ++i) {
        res[i] -= p2[i];
    }
    return res;
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::operator*(const Point<N, T>& point, T scalar) {
    Point<N, T> res;
    for (size_t i = 0; i < N; ++i) {
        res[i] = point[i] * scalar;
    }
    return res;
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::operator*(T scalar, const Point<N, T>& point) {
    return point * scalar;
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::operator/(const Point<N, T>& point, T scalar) {
    return point * (1. / scalar);
}
