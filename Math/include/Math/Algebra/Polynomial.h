///@file Barrier.h
#pragma once

#include "math_export.hpp"

#include <vector>
#include <map>

namespace math::algebra {

    /**
        @brief Класс для работы с полиномами одной переменной.
        @details Представляет полином вида:
            P(x) = aₙxⁿ + aₙ₋₁xⁿ⁻¹ + ... + a₁x + a₀
            где aᵢ - коэффициенты при соответствующих степенях x.

        @note Поддерживает основные алгебраические операции:
            - Сложение, вычитание и умножение полиномов
            - Вычисление производной
            - Деление полиномов с остатком
            - Вычисление НОД двух полиномов
            - Вычисление значения в точке

        @note Коэффициенты хранятся в виде отображения "степень → коэффициент",
            что позволяет эффективно работать с разреженными полиномами.

        @warning При операциях с полиномами высокой степени возможно переполнение.
    */
    class MATH_EXPORT Polynomial {
    public:
        using value_type = double;  ///< Тип коэффициентов полинома

    public:
        Polynomial() = default;     ///< Создает нулевой полином (P(x) = 0)
        Polynomial(const std::vector<value_type>& coeffs);      ///<  Создает полином из вектора коэффициентов, где coeffs[i] соответствует коэффициенту при xⁱ
        Polynomial(const std::map<int, value_type>& coeffs);    ///< Создает полином из отображения "степень -> коэффициент". Пример: {{0,1}, {2,3}} создаст полином 3x² + 1
        virtual ~Polynomial() = default;    

        bool operator==(const Polynomial& other) const;     
        bool operator!=(const Polynomial& other) const;     

        Polynomial& operator+=(const Polynomial& other);    
        Polynomial& operator-=(const Polynomial& other);    
        Polynomial& operator*=(const Polynomial& other);    
        Polynomial derivative() const;      ///< Вычисляет производную полинома
        std::pair<Polynomial, Polynomial> divide(const Polynomial& other) const;    ///< Делит полином на другой с остатком

        value_type evaluate(value_type x) const;    ///< Вычисляет значение полинома в точке x
        value_type get_coeff(int degree) const;     
        void set_coeff(int degree, value_type value);
        int degree() const;     ///< Возвращает степень полинома 

        static Polynomial gcd(const Polynomial& p1, const Polynomial& p2);  ///< Вычисляет наибольший общий делитель двух полиномов, используя алгоритм Евклида

    protected:
        std::map<int, value_type> m_coeffs;     ///< Хранилище коэффициентов (степень -> коэффициент)
    };

    MATH_EXPORT Polynomial operator+(const Polynomial& p1, const Polynomial& p2);
    MATH_EXPORT Polynomial operator-(const Polynomial& p1, const Polynomial& p2);
    MATH_EXPORT Polynomial operator*(const Polynomial& p1, const Polynomial& p2);

} // namespace math::algebra
