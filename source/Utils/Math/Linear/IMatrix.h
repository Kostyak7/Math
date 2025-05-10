#pragma once 

#include "DVector.h"

#include <complex>

namespace math::linal {

    class IMatrix {
    public:
        using value_type = double;
        using complex_value_type = std::complex<value_type>;

    public:
        virtual ~IMatrix() = default;

        virtual bool is_equal(const IMatrix& other) const;
        virtual bool is_empty() const;

        virtual bool is_square() const;
        virtual bool is_zero() const;
        virtual bool is_identity() const;
        virtual bool is_diagonal() const;
        virtual bool is_symmetrical() const;
        virtual bool is_upper_triangular() const;
        virtual bool is_lower_triangular() const;
        virtual bool is_trapezoidal() const;
        virtual bool is_positive_definite() const;
        virtual bool is_negative_definite() const;

        virtual size_t get_width() const = 0;
        virtual size_t get_height() const = 0;

        virtual value_type get(size_t row, size_t col) const = 0;
        virtual void set(size_t row, size_t col, value_type value) = 0;

        virtual value_type det() const = 0;

        virtual std::vector<std::pair<size_t, complex_value_type>> get_eigenvalues() const = 0;
        virtual std::vector<std::pair<complex_value_type, std::vector<DVector>>> get_eigenvectors() const = 0;
    };

    bool is_positive_definite_stochastic(const IMatrix& matrix, size_t test_points = 1000);
    bool is_negative_definite_stochastic(const IMatrix& matrix, size_t test_points = 1000);

} // namespace math::linal