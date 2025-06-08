#pragma once

#include "math_export.hpp"

namespace math::geom {

    class MATH_EXPORT IVector {
    public:
        using value_type = double;

        const value_type& front() const noexcept;
        value_type& front() noexcept;
        const value_type& back() const noexcept;
        value_type& back() noexcept;

    public:
        ~IVector() = default;

        virtual size_t size() const = 0;
        virtual bool empty() const;

        void fill(value_type value);

        bool is_zero() const;

        const value_type& operator[](size_t pos) const noexcept;
        value_type& operator[](size_t pos) noexcept;

        virtual const value_type& at(size_t pos) const = 0;
        virtual value_type& at(size_t pos) = 0;

        value_type x() const;
        value_type y() const;
        value_type z() const;
        value_type w() const;

        value_type dot(const IVector& other) const;
        value_type norm() const;
        value_type length() const;

    protected:
        virtual const value_type& _(size_t pos) const noexcept = 0;
        virtual value_type& _(size_t pos) noexcept = 0;
    };

} // namespace math::geom
