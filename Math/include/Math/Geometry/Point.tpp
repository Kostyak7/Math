template <size_t N, typename T>
typename math::geom::Point<N, T>::iterator math::geom::Point<N, T>::begin() noexcept {
    return m_data.begin();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::const_iterator math::geom::Point<N, T>::begin() const noexcept {
    return m_data.begin();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::iterator math::geom::Point<N, T>::end() noexcept {
    return m_data.end();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::const_iterator math::geom::Point<N, T>::end() const noexcept {
    return m_data.end();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::reverse_iterator math::geom::Point<N, T>::rbegin() noexcept {
    return m_data.rbegin();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::const_reverse_iterator math::geom::Point<N, T>::rbegin() const noexcept {
    return m_data.rbegin();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::reverse_iterator math::geom::Point<N, T>::rend() noexcept {
    return m_data.rend();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::const_reverse_iterator math::geom::Point<N, T>::rend() const noexcept {
    return m_data.rend();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::const_iterator math::geom::Point<N, T>::cbegin() const noexcept {
    return m_data.cbegin();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::const_iterator math::geom::Point<N, T>::cend() const noexcept {
    return m_data.cend();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::const_reverse_iterator math::geom::Point<N, T>::crbegin() const noexcept {
    return m_data.crbegin();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::const_reverse_iterator math::geom::Point<N, T>::crend() const noexcept {
    return m_data.crend();
}

template <size_t N, typename T>
typename const math::geom::Point<N, T>::value_type& math::geom::Point<N, T>::front() const noexcept {
    return m_data.front();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type& math::geom::Point<N, T>::front() noexcept {
    return m_data.front();
}

template <size_t N, typename T>
typename const math::geom::Point<N, T>::value_type& math::geom::Point<N, T>::back() const noexcept {
    return m_data.back();
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type& math::geom::Point<N, T>::back() noexcept {
    return m_data.back();
}

template <size_t N, typename T>
math::geom::Point<N, T>::Point(value_type default_value) {
    std::fill(m_data.begin(), m_data.end(), default_value);
}

template <size_t N, typename T>
math::geom::Point<N, T>::Point(std::initializer_list<value_type> init) {
    if (init.size() != N) {
        // throw ?
        return;
    }
    std::copy(init.begin(), init.end(), begin());
}

template <size_t N, typename T>
math::geom::Point<N, T>::Point(const base_container_type& container) {
    if (container.size() != N) {
        // throw ?
        return;
    }
    std::copy(container.begin(), container.end(), begin());
}

template <size_t N, typename T>
math::geom::Point<N, T>::Point(const Point& other)
    : m_data(other.m_data)
{
}

template <size_t N, typename T>
template <class T2>
math::geom::Point<N, T>::Point(const Point<N, T2>& other) {
    for (size_t i = 0; i < N; ++i) {
        m_data[i] = static_cast<value_type>(other.m_data[i]);
    }
}

template <size_t N, typename T>
math::geom::Point<N, T>::Point(Point&& other) noexcept
    : m_data(std::move(m_data))
{
}

template <size_t N, typename T>
math::geom::Point<N, T>& math::geom::Point<N, T>::operator=(const Point& other) {
    if (this == &other) return *this;
    for (size_t i = 0; i < N; ++i) {
        m_data[i] = other.m_data[i];
    }
    return *this;
}

template <size_t N, typename T>
template <class T2>
math::geom::Point<N, T>& math::geom::Point<N, T>::operator=(const Point<N, T2>& other) {
    if (this == &other) return *this;
    for (size_t i = 0; i < N; ++i) {
        m_data[i] = static_cast<T2>(other.m_data[i]);
    }
    return *this;
}

template <size_t N, typename T>
math::geom::Point<N, T>& math::geom::Point<N, T>::operator=(Point&& other) noexcept {
    if (this == &other) return *this;
    m_data = std::move(other.m_data);
    return *this;
}

template <size_t N, typename T>
math::geom::Point<N, T>& math::geom::Point<N, T>::operator*=(value_type scalar) {
    for (size_t i = 0; i < N; ++i) {
        m_data[i] *= scalar;
    }
    return *this;
}

template <size_t N, typename T>
math::geom::Point<N, T>& math::geom::Point<N, T>::operator/=(value_type scalar) {
    for (size_t i = 0; i < N; ++i) {
        m_data[i] /= scalar;
    }
    return *this;
}

template <size_t N, typename T>
math::geom::Point<N, T>& math::geom::Point<N, T>::operator+=(const Point& other) {
    for (size_t i = 0; i < N; ++i) {
        m_data[i] += other.m_data[i];
    }
    return *this;
}

template <size_t N, typename T>
math::geom::Point<N, T>& math::geom::Point<N, T>::operator-=(const Point& other) {
    for (size_t i = 0; i < N; ++i) {
        m_data[i] -= other.m_data[i];
    }
    return *this;
}

template <size_t N, typename T>
void math::geom::Point<N, T>::swap(Point& other) noexcept {
    std::swap(m_data, other.m_data);
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type math::geom::Point<N, T>::distance_to(const Point& other) {
    value_type res2 = 0;
    for (size_t i = 0; i < N; ++i) {
        const auto di = m_data[i] - other.m_data[i];
        res2 += di * di;
    }
    return std::sqrt(res2);
}

template <size_t N, typename T>
size_t math::geom::Point<N, T>::size() const {
    return N;
}

template <size_t N, typename T>
bool math::geom::Point<N, T>::empty() const {
    return N == 0;
}

template <size_t N, typename T>
void math::geom::Point<N, T>::fill(value_type value) {
    m_data.fill(value);
}

template <size_t N, typename T>
bool math::geom::Point<N, T>::is_zero() const {
    for (size_t i = 0; i < N; ++i) {
        if (dcmp(m_data[i]) != 0) {
            return false;
        }
    }
    return true;
}

template <size_t N, typename T>
typename const math::geom::Point<N, T>::value_type& math::geom::Point<N, T>::operator[](size_t pos) const noexcept {
    return m_data[pos];
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type& math::geom::Point<N, T>::operator[](size_t pos) noexcept {
    return m_data[pos];
}

template <size_t N, typename T>
typename const math::geom::Point<N, T>::value_type& math::geom::Point<N, T>::at(size_t pos) const {
    return m_data.at(pos);
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type& math::geom::Point<N, T>::at(size_t pos) {
    return m_data.at(pos);
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type math::geom::Point<N, T>::x() const {
    if constexpr (N >= 1) {
        return m_data[0];
    }
    return {};
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type math::geom::Point<N, T>::y() const {
    if constexpr (N >= 2) {
        return m_data[1];
    }
    return {};
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type math::geom::Point<N, T>::z() const {
    if constexpr (N >= 3) {
        return m_data[2];
    }
    return {};
}

template <size_t N, typename T>
typename math::geom::Point<N, T>::value_type math::geom::Point<N, T>::w() const {
    if constexpr (N >= 4) {
        return m_data[3];
    }
    return {};
}

template <size_t N, typename T>
void math::geom::swap(Point<N, T>& p1, Point<N, T>& p2) noexcept {
    p1.swap(p2);
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::operator*(Point<N, T> point, T scalar) {
    Point<N, T> res;
    for (size_t i = 0; i < N; ++i) {
        res[i] = point[i] * scalar;
    }
    return res;
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::operator*(T scalar, Point<N, T> point) {
    return point * scalar;
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::operator/(Point<N, T> point, T scalar) {
    return point * (1. / scalar);
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::operator+(Point<N, T> p1, const Point<N, T>& p2) {
    for (size_t i = 0; i < N; ++i) {
        p1[i] += p2[i];
    }
    return p1;
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::operator-(Point<N, T> p1, const Point<N, T>& p2) {
    for (size_t i = 0; i < N; ++i) {
        p1[i] -= p2[i];
    }
    return p1;
}

template <size_t N, typename T>
typename math::geom::Point<N, T> math::geom::operator-(const Point<N, T>& point) {
    Point<N, T> res;
    for (size_t i = 0; i < N; ++i) {
        res[i] = -point[i];
    }
    return res;
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::zero_point() {
    return {};
}

template <size_t N, typename T>
math::geom::Point<N, T> math::geom::unit_point(size_t i, T value) {
    Point<N, T> res;
    if (i < N) {
        res[i] = value;
    }
    return res;
}
