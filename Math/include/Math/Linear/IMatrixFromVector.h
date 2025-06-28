#pragma once 

#include "IMatrix.h"

namespace math::linal {

    class MATH_EXPORT IMatrixFromVector : public IMatrix {
    public:
        class ConstProxyVector;
        class ProxyVector;

    public:
        IMatrixFromVector() noexcept = default;
        IMatrixFromVector(size_t size, const value_type& default_value = {});
        IMatrixFromVector(std::vector<value_type> data);

        ConstProxyVector operator[](size_t row) const;
        ProxyVector operator[](size_t row);

        ConstProxyVector at(size_t row) const;
        ProxyVector at(size_t row);

        void clear() override;

    protected:
        virtual ConstProxyVector get_row(size_t row) const = 0;
        virtual ProxyVector get_row(size_t row) = 0;
        virtual bool check_row_index(size_t row) const = 0;
        virtual size_t calculate_index(size_t row, size_t col) const = 0;
        
        const value_type& _(size_t row, size_t col) const; //direct access
        value_type& _(size_t row, size_t col); // Its not an operator() because in methods it is perceived as an operator,

    protected:
        std::vector<value_type> m_data;
    };

    class IMatrixFromVector::ConstProxyVector {
    public:
        using const_iterator = typename std::vector<value_type>::const_iterator;

        const_iterator begin() const;
        const_iterator end() const;

    public:
        ConstProxyVector(size_t size, const_iterator front);

        value_type operator[](size_t col) const;
        value_type at(size_t col) const;

        size_t size() const;
        bool empty() const;

        operator DVector() const;

    private:
        size_t m_size;
        const_iterator m_front;
    };

    MATH_EXPORT bool operator==(const IMatrixFromVector::ConstProxyVector& v1, const IMatrixFromVector::ConstProxyVector& v2);
    MATH_EXPORT bool operator!=(const IMatrixFromVector::ConstProxyVector& v1, const IMatrixFromVector::ConstProxyVector& v2);

    class IMatrixFromVector::ProxyVector {
    public:
        using const_iterator = typename std::vector<value_type>::const_iterator;
        using iterator = typename std::vector<value_type>::iterator;

        const_iterator begin() const;
        const_iterator end() const;
        iterator begin();
        iterator end();

    public:
        ProxyVector(size_t size, iterator front);

        value_type operator[](size_t col) const;
        value_type& operator[](size_t col);

        value_type at(size_t col) const;
        value_type& at(size_t col);

        size_t size() const;
        bool empty() const;

        operator DVector() const;

    private:
        size_t m_size;
        iterator m_front;
    };

    MATH_EXPORT bool operator==(const IMatrixFromVector::ProxyVector& v1, const IMatrixFromVector::ProxyVector& v2);
    MATH_EXPORT bool operator!=(const IMatrixFromVector::ProxyVector& v1, const IMatrixFromVector::ProxyVector& v2);

} // namespace math::linal
