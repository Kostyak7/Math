template <size_t N>
math::geom::Vector<N>::Vector(const Point<N, value_type>& to) {
    std::copy(to.begin(), to.end(), this->begin());
}

template <size_t N>
math::geom::Vector<N>::Vector(const Point<N, value_type>& from, const Point<N, value_type>& to)
    : Vector(to - from)
{
}

template <size_t N>
math::geom::Vector<N>::Vector(std::initializer_list<value_type> init) {
    std::copy(init.begin(), init.end(), this->begin());
}

template <size_t N>
math::geom::Vector<N>& math::geom::Vector<N>::operator*=(value_type scalar) {
    for (auto& val : *this) {
        val *= scalar;
    }
    return *this;
}

template <size_t N>
math::geom::Vector<N>& math::geom::Vector<N>::operator/=(value_type scalar) {
    return *this *= (1. / scalar);
}

template <size_t N>
math::geom::Vector<N>& math::geom::Vector<N>::operator+=(const Vector& other) {
    for (size_t i = 0; i < N; ++i) {
        (*this)[i] += other[i];
    }
    return *this;
}

template <size_t N>
math::geom::Vector<N>& math::geom::Vector<N>::operator-=(const Vector& other) {
    for (size_t i = 0; i < N; ++i) {
        (*this)[i] -= other[i];
    }
    return *this;
}

template <size_t N>
bool math::geom::Vector<N>::is_zero() const {
    for (const auto& val : *this) {
        if (dcmp(val) != 0) {
            return false;
        }
    }
    return true;
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::x() const {
    if constexpr (N >= 1) {
        return (*this)[0];
    }
    return {};
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::y() const {
    if constexpr (N >= 2) {
        return (*this)[1];
    }
    return {};
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::z() const {
    if constexpr (N >= 3) {
        return (*this)[2];
    }
    return {};
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::w() const {
    if constexpr (N >= 4) {
        return (*this)[3];
    }
    return {};
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::dot(const Vector& other) const {
    value_type res = 0;
    for (size_t i = 0; i < N; ++i) {
        res += (*this)[i] * other[i];
    }
    return res;
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::norm() const {
    return std::sqrt(dot(*this));
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::length() const {
    return norm();
}

template <size_t N>
math::geom::Vector<N>& math::geom::Vector<N>::normalize() {
    if (value_type len = norm(); dcmp(len) != 0) {
        for (size_t i = 0; i < N; ++i) {
            (*this)[i] /= len;
        }
    }
    return *this;
}

template <size_t N>
math::geom::Vector<N> math::geom::operator*(const Vector<N>& vec, typename  Vector<N>::value_type scalar) {
    Vector<N> res(vec);
    res *= scalar;
    return res;
}

template <size_t N>
math::geom::Vector<N> math::geom::operator*(typename Vector<N>::value_type scalar, const Vector<N>& vec) {
    return vec * scalar;
}

template <size_t N>
math::geom::Vector<N> math::geom::operator/(const Vector<N>& vec, typename Vector<N>::value_type scalar) {
    return vec * (1. / scalar);
}

template <size_t N>
math::geom::Vector<N> math::geom::operator+(const Vector<N>& v1, const Vector<N>& v2) {
    Vector<N> result(v1);
    return result += v2;
}

template <size_t N>
math::geom::Vector<N> math::geom::operator-(const Vector<N>& v1, const Vector<N>& v2) {
    Vector<N> result(v1);
    return result -= v2;
}

template <size_t N>
math::geom::Vector<N> math::geom::operator-(const Vector<N>& vec) {
    Vector<N> res;
    for (size_t i = 0; i < N; ++i) {
        res[i] = -vec[i];
    }
    return res;
}

template <size_t N>
math::geom::Vector<N> math::geom::normalized(const Vector<N>& vec) {
    Vector<N> res(vec);
    res.normalize();
    return res;
}

typename math::geom::Vector2D::value_type math::geom::pseudodot(const Vector2D& v1, const Vector2D& v2) {
    return v1.x() * v2.y() - v1.y() * v2.x();
}

typename math::geom::Vector3D math::geom::cross_product(const Vector3D& v1, const Vector3D& v2) {
    return {
        v1.y() * v2.z() - v1.z() * v2.y(),
        v1.z() * v2.x() - v1.x() * v2.z(),
        v1.x() * v2.y() - v1.y() * v2.x(),
    };
}
