template <size_t N>
typename math::geom::Vector<N>::iterator math::geom::Vector<N>::begin() noexcept {
    return m_data.begin();
}

template <size_t N>
typename math::geom::Vector<N>::const_iterator math::geom::Vector<N>::begin() const noexcept {
    return m_data.begin();
}

template <size_t N>
typename math::geom::Vector<N>::iterator math::geom::Vector<N>::end() noexcept {
    return m_data.end();
}

template <size_t N>
typename math::geom::Vector<N>::const_iterator math::geom::Vector<N>::end() const noexcept {
    return m_data.end();
}

template <size_t N>
typename math::geom::Vector<N>::reverse_iterator math::geom::Vector<N>::rbegin() noexcept {
    return m_data.rbegin();
}

template <size_t N>
typename math::geom::Vector<N>::const_reverse_iterator math::geom::Vector<N>::rbegin() const noexcept {
    return m_data.rbegin();
}

template <size_t N>
typename math::geom::Vector<N>::reverse_iterator math::geom::Vector<N>::rend() noexcept {
    return m_data.rend();
}

template <size_t N>
typename math::geom::Vector<N>::const_reverse_iterator math::geom::Vector<N>::rend() const noexcept {
    return m_data.rend();
}

template <size_t N>
typename math::geom::Vector<N>::const_iterator math::geom::Vector<N>::cbegin() const noexcept {
    return m_data.cbegin();
}

template <size_t N>
typename math::geom::Vector<N>::const_iterator math::geom::Vector<N>::cend() const noexcept {
    return m_data.cend();
}

template <size_t N>
typename math::geom::Vector<N>::const_reverse_iterator math::geom::Vector<N>::crbegin() const noexcept {
    return m_data.crbegin();
}

template <size_t N>
typename math::geom::Vector<N>::const_reverse_iterator math::geom::Vector<N>::crend() const noexcept {
    return m_data.crend();
}

template <size_t N>
typename const math::geom::Vector<N>::value_type& math::geom::Vector<N>::front() const noexcept {
    return m_data.front();
}

template <size_t N>
typename math::geom::Vector<N>::value_type& math::geom::Vector<N>::front() noexcept {
    return m_data.front();
}

template <size_t N>
typename const math::geom::Vector<N>::value_type& math::geom::Vector<N>::back() const noexcept {
    return m_data.back();
}

template <size_t N>
typename math::geom::Vector<N>::value_type& math::geom::Vector<N>::back() noexcept {
    return m_data.back();
}

template <size_t N>
math::geom::Vector<N>::Vector(value_type default_value) {
    std::fill(m_data.begin(), m_data.end(), default_value);
}

template <size_t N>
math::geom::Vector<N>::Vector(const Point<N, value_type>& to) {
    std::copy(to.begin(), to.end(), begin());
}

template <size_t N>
math::geom::Vector<N>::Vector(const Point<N, value_type>& from, const Point<N, value_type>& to)
    : Vector(to - from)
{
}

template <size_t N>
math::geom::Vector<N>::Vector(std::initializer_list<value_type> init) {
    if (init.size() != N) {
        // throw ?
        return;
    }
    std::copy(init.begin(), init.end(), begin());
}

template <size_t N>
math::geom::Vector<N>::Vector(const base_container_type& container) {
    if (container.size() != N) {
        // throw ?
        return;
    }
    std::copy(container.begin(), container.end(), begin());
}

template <size_t N>
math::geom::Vector<N>::Vector(const Vector& other)
    : m_data(other.m_data)
{
}

template <size_t N>
math::geom::Vector<N>::Vector(Vector&& other) noexcept
    : m_data(std::move(other.m_data))
{
}

template <size_t N>
math::geom::Vector<N>& math::geom::Vector<N>::operator=(const Vector& other) {
    if (this == &other) return *this;
    m_data = other.m_data;
    return *this;
}

template <size_t N>
math::geom::Vector<N>& math::geom::Vector<N>::operator=(Vector&& other) noexcept {
    if (this == &other) return *this;
    m_data = std::move(other.m_data);
    return *this;
}

template <size_t N>
math::geom::Vector<N>& math::geom::Vector<N>::operator*=(value_type scalar) {
    for (auto& val : m_data) {
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
        m_data[i] += other[i];
    }
    return *this;
}

template <size_t N>
math::geom::Vector<N>& math::geom::Vector<N>::operator-=(const Vector& other) {
    for (size_t i = 0; i < N; ++i) {
        m_data[i] -= other[i];
    }
    return *this;
}

template <size_t N>
void math::geom::Vector<N>::swap(Vector& other) noexcept {
    std::swap(m_data, other.m_data);
}

template <size_t N>
math::geom::Point<N> math::geom::Vector<N>::get_end_point() const {
    return m_data;
}

template <size_t N>
size_t math::geom::Vector<N>::size() const {
    return N;
}

template <size_t N>
bool math::geom::Vector<N>::empty() const {
    return N == 0;
}

template <size_t N>
bool math::geom::Vector<N>::is_zero() const {
    for (const auto& val : m_data) {
        if (dcmp(val) != 0) {
            return false;
        }
    }
    m_data.back();
    return true;
}

template <size_t N>
typename const math::geom::Vector<N>::value_type& math::geom::Vector<N>::operator[](size_t pos) const noexcept {
    return m_data[pos];
}

template <size_t N>
typename math::geom::Vector<N>::value_type& math::geom::Vector<N>::operator[](size_t pos) noexcept {
    return m_data[pos];
}

template <size_t N>
typename const math::geom::Vector<N>::value_type& math::geom::Vector<N>::at(size_t pos) const {
    return m_data.at(pos);
}

template <size_t N>
typename math::geom::Vector<N>::value_type& math::geom::Vector<N>::at(size_t pos) {
    return m_data.at(pos);
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::x() const {
    if constexpr (N >= 1) {
        return m_data[0];
    }
    return {};
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::y() const {
    if constexpr (N >= 2) {
        return m_data[1];
    }
    return {};
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::z() const {
    if constexpr (N >= 3) {
        return m_data[2];
    }
    return {};
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::w() const {
    if constexpr (N >= 4) {
        return m_data[3];
    }
    return {};
}

template <size_t N>
typename math::geom::Vector<N>::value_type math::geom::Vector<N>::dot(const Vector& other) const {
    value_type res = 0;
    for (size_t i = 0; i < N; ++i) {
        res += m_data[i] * other[i];
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
            m_data[i] /= len;
        }
    }
    return *this;
}

template <size_t N>
void math::geom::swap(Vector<N>& v1, Vector<N>& v2) noexcept {
    v1.swap(v2);
}

template <size_t N>
math::geom::Vector<N> math::geom::operator*(Vector<N> vec, typename  Vector<N>::value_type scalar) {
    vec *= scalar;
    return vec;
}

template <size_t N>
math::geom::Vector<N> math::geom::operator*(typename Vector<N>::value_type scalar, Vector<N> vec) {
    return vec * scalar;
}

template <size_t N>
math::geom::Vector<N> math::geom::operator/(Vector<N> vec, typename Vector<N>::value_type scalar) {
    return vec * (1. / scalar);
}

template <size_t N>
math::geom::Vector<N> math::geom::operator+(Vector<N> v1, const Vector<N>& v2) {
    v1 += v2;
    return v1;
}

template <size_t N>
math::geom::Vector<N> math::geom::operator-(Vector<N> v1, const Vector<N>& v2) {
    v1 -= v2;
    return v1;
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
math::geom::Vector<N> math::geom::normalized(Vector<N> vec) {
    vec.normalize();
    return vec;
}

template <size_t N>
math::geom::Vector<N> math::geom::zero_vector() {
    return {};
}

template <size_t N>
math::geom::Vector<N> math::geom::unit_vector(size_t i, typename Vector<N>::value_type value) {
    Vector<N> res{};
    if (i < N) {
        res[i] = value;
    }
    return res;
}

typename math::geom::Vector2D::value_type math::geom::pseudodot(const Vector2D& v1, const Vector2D& v2) {
    return v1.x() * v2.y() - v1.y() * v2.x();
}

math::geom::Vector3D math::geom::cross_product(const Vector3D& v1, const Vector3D& v2) {
    return {
        v1.y() * v2.z() - v1.z() * v2.y(),
        v1.z() * v2.x() - v1.x() * v2.z(),
        v1.x() * v2.y() - v1.y() * v2.x(),
    };
}
