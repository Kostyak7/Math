#include <Utils/Math/Algebra/Polynomial.h>

#include <Utils/Math/Common.h>

math::algebra::Polynomial::Polynomial(const std::vector<value_type>& coeffs) {
    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (dcmp(coeffs[i]) != 0) {
            m_coeffs.insert({ i, coeffs[i] });
        }
    }
}

math::algebra::Polynomial::Polynomial(const std::map<int, value_type>& coeffs)
    : m_coeffs(coeffs)
{
}

bool math::algebra::Polynomial::operator==(const Polynomial& other) const {
    if (degree() != other.degree()) {
        return false;
    }

    for (const auto& [degree, value] : m_coeffs) {
        auto it = other.m_coeffs.find(degree);
        if (it == other.m_coeffs.end()) {
            if (dcmp(value) != 0) {
                return false;
            }
        }
        else if (dcmp(value, it->second) != 0) {
            return false;
        }
    }

    for (const auto& [degree, value] : other.m_coeffs) {
        if (auto it = m_coeffs.find(degree); it == m_coeffs.end() && dcmp(value) != 0) {
            return false;
        }
    }

    return true;
}

bool math::algebra::Polynomial::operator!=(const Polynomial& other) const {
    return !(*this == other);
}

math::algebra::Polynomial& math::algebra::Polynomial::operator+=(const Polynomial& other) {
    for (const auto& [deg, coef] : other.m_coeffs) {
        m_coeffs[deg] += coef;
        if (dcmp(m_coeffs[deg]) == 0)
            m_coeffs.erase(deg);
    }
    return *this;
}

math::algebra::Polynomial& math::algebra::Polynomial::operator-=(const Polynomial& other) {
    for (const auto& [deg, coef] : other.m_coeffs) {
        m_coeffs[deg] -= coef;
        if (dcmp(m_coeffs[deg]) == 0)
            m_coeffs.erase(deg);
    }
    return *this;
}

math::algebra::Polynomial& math::algebra::Polynomial::operator*=(const Polynomial& other) {
    if (other.m_coeffs.empty()) {
        m_coeffs.clear();
        return *this;
    }

    for (const auto& [deg1, coef1] : m_coeffs) {
        for (const auto& [deg2, coef2] : other.m_coeffs) {
            m_coeffs[deg1 + deg2] += coef1 * coef2;
        }
    }
    return *this;
}

math::algebra::Polynomial math::algebra::Polynomial::derivative() const {
    Polynomial res;
    for (const auto& [deg, coef] : m_coeffs) {
        if (deg > 0) {
            res.m_coeffs[deg - 1] = deg * coef;
        }
    }
    return res;
}

std::pair<math::algebra::Polynomial, math::algebra::Polynomial> math::algebra::Polynomial::divide(const Polynomial& other) const {
    if (other.m_coeffs.empty()) {
        throw std::invalid_argument("Division by zero polynomial");
    }

    Polynomial quotient;
    Polynomial remainder = *this;

    while (!remainder.m_coeffs.empty() && remainder.degree() >= other.degree()) {
        int deg_diff = remainder.degree() - other.degree();
        value_type coef_ratio = remainder.get_coeff(remainder.degree()) / other.get_coeff(other.degree());

        Polynomial term(std::map<int, value_type>{ {deg_diff, coef_ratio} });
        quotient = quotient + term;
        remainder = remainder - (term * other);
    }

    return { quotient, remainder };
}

math::algebra::Polynomial::value_type math::algebra::Polynomial::evaluate(value_type x) const {
    value_type result = 0;
    for (const auto& [deg, coef] : m_coeffs) {
        result += coef * std::pow(x, deg);
    }
    return result;
}

math::algebra::Polynomial::value_type math::algebra::Polynomial::get_coeff(int degree) const {
    if (auto it = m_coeffs.find(degree); it != m_coeffs.end()) {
        return it->second;
    }
    return 0.0;
}

void math::algebra::Polynomial::set_coeff(int degree, value_type value) {
    if (auto it = m_coeffs.find(degree); it != m_coeffs.end()) {
        it->second += value;
    }
    else {
        m_coeffs.insert({ degree, value });
    }

}

int math::algebra::Polynomial::degree() const {
    if (m_coeffs.empty()) return -1;
    return m_coeffs.rbegin()->first;
}

math::algebra::Polynomial math::algebra::Polynomial::gcd(const Polynomial& p1, const Polynomial& p2) {
    if (p2.m_coeffs.empty()) return p1;
    return gcd(p2, p1.divide(p2).second);
}

math::algebra::Polynomial math::algebra::operator+(const Polynomial& p1, const Polynomial& p2) {
    Polynomial res = p1;
    res += p2;
    return res;
}

math::algebra::Polynomial math::algebra::operator-(const Polynomial& p1, const Polynomial& p2) {
    Polynomial res = p1;
    res -= p2;
    return res;
}

math::algebra::Polynomial math::algebra::operator*(const Polynomial& p1, const Polynomial& p2) {
    Polynomial res = p1;
    res *= p2;
    return res;
}
