#pragma once 

#include "math_export.hpp"

#include <vector>

namespace math::linal {

	class MATH_EXPORT DVector {
    public:
        using value_type = double;
        using base_container_type = std::vector<value_type>;
        using const_iterator = typename base_container_type::const_iterator;
        using iterator = typename base_container_type::iterator;
        using const_reverse_iterator = typename base_container_type::const_reverse_iterator;
        using reverse_iterator = typename base_container_type::reverse_iterator;

        iterator begin() noexcept;
        const_iterator begin() const noexcept;
        iterator end() noexcept;
        const_iterator end() const noexcept;
        reverse_iterator rbegin() noexcept;
        const_reverse_iterator rbegin() const noexcept;
        reverse_iterator rend() noexcept;
        const_reverse_iterator rend() const noexcept;
        const_iterator cbegin() const noexcept;
        const_iterator cend() const noexcept;
        const_reverse_iterator crbegin() const noexcept;
        const_reverse_iterator crend() const noexcept;

        const value_type& front() const noexcept;
        value_type& front() noexcept;
        const value_type& back() const noexcept;
        value_type& back() noexcept;

	public:
        DVector() = default;
        DVector(size_t size, value_type default_val = 0.0);
        DVector(std::initializer_list<value_type> init);
        DVector(const base_container_type& container);
        DVector(const DVector& other);
        DVector(DVector&& other) noexcept;
        ~DVector() = default;

        DVector& operator=(const DVector& other);
        DVector& operator=(DVector&& other) noexcept;

        DVector& operator*=(value_type scalar);
        DVector& operator/=(value_type scalar);

        DVector& operator+=(const DVector& other);
        DVector& operator-=(const DVector& other);

        void swap(DVector& other) noexcept;

        size_t size() const;
        bool empty() const;
        void resize(size_t new_size);
        void resize(size_t new_size, value_type val);
        void fill(value_type value);

        bool is_zero() const;

        const value_type& operator[](size_t pos) const noexcept;
        value_type& operator[](size_t pos) noexcept;

        const value_type& at(size_t pos) const;
        value_type& at(size_t pos);

        value_type dot(const DVector& other) const;
        value_type norm() const;
        value_type length() const;

        DVector& normalize();

    private:
        base_container_type m_data;
	};

    void MATH_EXPORT swap(DVector& v1, DVector& v2) noexcept;

    MATH_EXPORT bool operator==(const DVector& v1, const DVector& v2);
    MATH_EXPORT bool operator!=(const DVector& v1, const DVector& v2);

    MATH_EXPORT DVector operator*(DVector vector, DVector::value_type scalar);
    MATH_EXPORT DVector operator*(DVector::value_type scalar, DVector vector);
    MATH_EXPORT DVector operator/(DVector vector, DVector::value_type scalar);

    MATH_EXPORT DVector operator+(const DVector& v1, const DVector& v2);
    MATH_EXPORT DVector operator-(const DVector& v1, const DVector& v2);

    MATH_EXPORT DVector operator-(const DVector& vector);

    DVector MATH_EXPORT normalized(DVector vector);

} // namespace math::linal
