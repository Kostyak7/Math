#pragma once 

#include "IMatrix.h"

namespace math::linal {

    /**
     @brief Ленточная матрица особого вида
     @details Симметричная, хранятся только строки начиная с диагонального элемента
    */
    class BandMatrix final: public IMatrix, 
                            public std::vector<DVector> {
    public:
        using value_type = IMatrix::value_type;
        using std::vector<DVector>::vector;

    public:
        BandMatrix& operator*=(value_type scalar);
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

        size_t get_width() const override;
        size_t get_height() const override;

        size_t get_isl() const;
        void set_isl(size_t width);        
        void shrink();

        value_type det() const override;

        std::vector<std::pair<size_t, complex_value_type>> get_eigenvalues() const override;
        std::vector<std::pair<complex_value_type, std::vector<DVector>>> get_eigenvectors() const override;

        static BandMatrix identity_matrix(size_t n, size_t isl = 1);
        static BandMatrix elementary_matrix_unit(size_t n, size_t i, size_t j, double value = 1.0);
    };

    bool operator==(const BandMatrix& m1, const BandMatrix& m2);
    bool operator!=(const BandMatrix& m1, const BandMatrix& m2);

    BandMatrix operator*(const BandMatrix& matrix, BandMatrix::value_type scalar);
    BandMatrix operator*(BandMatrix::value_type scalar, const BandMatrix& matrix);
    BandMatrix operator/(const BandMatrix& matrix, BandMatrix::value_type scalar);

    DVector operator*(const BandMatrix& matrix, const DVector& vector);
    DVector operator*(const DVector& vector, const BandMatrix& matrix);

    BandMatrix operator*(const BandMatrix& m1, const BandMatrix& m2);

    BandMatrix operator+(const BandMatrix& m1, const BandMatrix& m2);
    BandMatrix operator-(const BandMatrix& m1, const BandMatrix& m2);

    BandMatrix operator-(const BandMatrix& matrix);

    BandMatrix inversed(const BandMatrix& matrix);

} // namespace math::linal
