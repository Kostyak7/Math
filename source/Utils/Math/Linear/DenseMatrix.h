#pragma once 

#include "IMatrix.h"

namespace math::linal {

    /**
     @brief Обычная плотная матрица 
    */
    class DenseMatrix final: public IMatrix {
    public:
        class ProxyVector {
        public:
            ProxyVector(FVector& vector);

            const value_type& operator[](size_t i) const;
            value_type& operator[](size_t i);

        private:
            FVector& m_vector;
        };

    public:
        DenseMatrix() noexcept = default;
        DenseMatrix(size_t height, size_t width, const value_type& default_value = {});
        DenseMatrix(size_t height, const FVector& default_value = {});
        DenseMatrix(const DenseMatrix& matrix);
        DenseMatrix(DenseMatrix&& matrix) noexcept;
        DenseMatrix(std::initializer_list<value_type> list);
        DenseMatrix(std::initializer_list<std::initializer_list<value_type>> list);

        DenseMatrix& operator=(const DenseMatrix& other);
        DenseMatrix& operator=(DenseMatrix&& other);

        DenseMatrix& operator*=(value_type scalar) noexcept;
        DenseMatrix& operator/=(value_type scalar);

        DenseMatrix& operator+=(const DenseMatrix& other);
        DenseMatrix& operator-=(const DenseMatrix& other);     

        bool is_positive_definite() const override;

        value_type get(size_t row, size_t col) const override;
        void set(size_t row, size_t col, value_type value) override;

        const FVector& operator[](size_t row) const;
        ProxyVector operator[](size_t row);

        size_t get_width() const override;
        size_t get_height() const override;

        void clear() override;

        value_type det() const override;

        std::vector<std::pair<size_t, complex_value_type>> get_eigenvalues() const override;
        std::vector<std::pair<complex_value_type, std::vector<FVector>>> get_eigenvectors() const override;

        void swap(DenseMatrix& other) noexcept;

        static DenseMatrix identity_matrix(size_t n);

        const std::vector<FVector>& data() const;
        std::vector<FVector>& data();

    private:
        std::vector<FVector> m_data;

    private:
        friend DenseMatrix inversed(const DenseMatrix& matrix);
    };

    void swap(DenseMatrix& m1, DenseMatrix& m2) noexcept;

    bool operator==(const DenseMatrix& m1, const DenseMatrix& m2);
    bool operator!=(const DenseMatrix& m1, const DenseMatrix& m2);

    DenseMatrix operator*(const DenseMatrix& matrix, DenseMatrix::value_type scalar);
    DenseMatrix operator*(DenseMatrix::value_type scalar, const DenseMatrix& matrix);
    DenseMatrix operator/(const DenseMatrix& matrix, DenseMatrix::value_type scalar);

    FVector operator*(const DenseMatrix& matrix, const FVector& vector);
    FVector operator*(const FVector& vector, const DenseMatrix& matrix);

    DenseMatrix operator*(const DenseMatrix& m1, const DenseMatrix& m2);

    DenseMatrix operator+(const DenseMatrix& m1, const DenseMatrix& m2);
    DenseMatrix operator-(const DenseMatrix& m1, const DenseMatrix& m2);

    DenseMatrix operator-(const DenseMatrix& matrix);

    DenseMatrix inversed(const DenseMatrix& matrix);

    DenseMatrix transposed(const DenseMatrix& matrix);

} // namespace math::linal