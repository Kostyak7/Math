#pragma once

#include <vector>
#include <map>

namespace math::algebra {

    class Polynomial {
    public:
        using value_type = double;

    public:
        Polynomial() = default;
        Polynomial(const std::vector<value_type>& coeffs);
        Polynomial(const std::map<int, value_type>& coeffs);
        virtual ~Polynomial() = default;

        bool operator==(const Polynomial& other) const;
        bool operator!=(const Polynomial& other) const;

        Polynomial& operator+=(const Polynomial& other) const;
        Polynomial& operator-=(const Polynomial& other) const;
        Polynomial& operator*=(const Polynomial& other) const;
        Polynomial derivative() const;
        std::pair<Polynomial, Polynomial> divide(const Polynomial& other) const;

        value_type evaluate(value_type x) const;
        value_type get_coeff(int degree) const;
        void set_coeff(int degree, value_type value);
        int degree() const;

        static Polynomial gcd(const Polynomial& p1, const Polynomial& p2);

    protected:
        std::map<int, value_type> m_coeffs;
    };

    Polynomial operator+(const Polynomial& p1, const Polynomial& p2);
    Polynomial operator-(const Polynomial& p1, const Polynomial& p2);
    Polynomial operator*(const Polynomial& p1, const Polynomial& p2);

} // namespace math::algebra
