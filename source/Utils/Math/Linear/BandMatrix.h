#pragma once 

#include "IMatrix.h"

namespace math::linal {

    /**
     @brief Ленточная матрица особого вида
     @details Симметричная, хранятся только строки начиная с диагонального элемента
    */
    class BandMatrix final: public IMatrix {
    public:
        BandMatrix() noexcept = default;
        BandMatrix(size_t n, const FVector& default_value = {});        
        BandMatrix(const BandMatrix& matrix);
        BandMatrix(BandMatrix&& matrix) noexcept;
        BandMatrix(std::initializer_list<std::initializer_list<value_type>> list);

        BandMatrix& operator=(const BandMatrix& other);
        BandMatrix& operator=(BandMatrix&& other);

        BandMatrix& operator*=(value_type scalar) noexcept;
        BandMatrix& operator/=(value_type scalar);

        BandMatrix& operator+=(const BandMatrix& other);
        BandMatrix& operator-=(const BandMatrix& other);

        bool is_zero() const override;
        bool is_identity() const override;
        bool is_diagonal() const override;
        bool is_symmetrical() const override;
        bool is_upper_triangular() const override;
        bool is_lower_triangular() const override;
        bool is_trapezoidal() const override;

        value_type get(size_t row, size_t col) const override;
        void set(size_t row, size_t col, value_type value) override;

        const FVector& operator[](size_t row) const;
        FVector& operator[](size_t row);

        size_t get_width() const override;
        size_t get_height() const override;

        size_t get_isl() const;
        void set_isl(size_t width);        
        void shrink();
        void clear() override;

        value_type det() const override;

        std::vector<std::pair<size_t, complex_value_type>> get_eigenvalues() const override;
        std::vector<std::pair<complex_value_type, std::vector<FVector>>> get_eigenvectors() const override;

        void swap(BandMatrix& other) noexcept;

        static BandMatrix identity_matrix(size_t n, size_t isl = 1);

        const std::vector<FVector>& data() const;
        std::vector<FVector>& data();

    private:
        std::vector<FVector> m_data;
    };

    void swap(BandMatrix& m1, BandMatrix& m2) noexcept;

    bool operator==(const BandMatrix& m1, const BandMatrix& m2);
    bool operator!=(const BandMatrix& m1, const BandMatrix& m2);

    BandMatrix operator*(const BandMatrix& matrix, BandMatrix::value_type scalar);
    BandMatrix operator*(BandMatrix::value_type scalar, const BandMatrix& matrix);
    BandMatrix operator/(const BandMatrix& matrix, BandMatrix::value_type scalar);

    FVector operator*(const BandMatrix& matrix, const FVector& vector);
    FVector operator*(const FVector& vector, const BandMatrix& matrix);

    BandMatrix operator*(const BandMatrix& m1, const BandMatrix& m2);

    BandMatrix operator+(const BandMatrix& m1, const BandMatrix& m2);
    BandMatrix operator-(const BandMatrix& m1, const BandMatrix& m2);

    BandMatrix operator-(const BandMatrix& matrix);

    BandMatrix inversed(const BandMatrix& matrix);

} // namespace math::linal