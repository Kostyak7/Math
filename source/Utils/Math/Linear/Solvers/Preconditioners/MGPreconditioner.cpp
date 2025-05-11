#include "MGPreconditioner.h"

#include <Utils/Math/Linear/Solvers/GaussSystemSolver.h>
#include <Utils/Math/Linear/Solvers/Preconditioners/ILUPreconditioner.h>
#include <Utils/Math/Common.h>

#include <algorithm>
#include <stdexcept>

class math::linal::MGPreconditioner::Impl : public math::linal::IPreconditioner::IImpl {
public:
    using MGType = MGPreconditioner::MGType;
    using Params = MGPreconditioner::Params;
    using SmootherType = MGPreconditioner::SmootherType;
    using CycleType = MGPreconditioner::CycleType;

    using Vector = math::linal::DVector;

public:
    Impl(const math::linal::AnyMatrix& matrix, const Params& params)
        : m_params(params)
        , m_matrix(matrix)
    {
    }

    Vector apply(const Vector& x) const override {
        if (check_levels()) return x;
        Vector result(x.size(), 0.0);
        switch (m_params.cycle_type) {
        case CycleType::V_CYCLE: v_cycle(0, x, result); break;
        case CycleType::W_CYCLE: w_cycle(0, x, result); break;
        case CycleType::F_CYCLE: f_cycle(0, x, result); break;
        }
        return result;
    }

protected:
    virtual bool check_levels() const = 0;

    virtual void v_cycle(size_t level_idx, const Vector& rhs, Vector& x) const = 0;
    virtual void w_cycle(size_t level_idx, const Vector& rhs, Vector& x) const = 0;
    virtual void f_cycle(size_t level_idx, const Vector& rhs, Vector& x) const = 0;

    virtual Vector restrict_algebraic(const Vector& fine, size_t level_idx) const = 0;
    virtual Vector restrict_geometric(const Vector& fine, size_t level_idx) const = 0;

    Vector restrict(const Vector& fine, size_t level_idx) const {
        switch (m_params.mg_type) {
                case MGType::Algebraic: return restrict_algebraic(fine, level_idx);
                case MGType::Geometric: return restrict_geometric(fine, level_idx);
            }
    }

    void coarse_sovler_init() {
        m_coarse_solver = std::make_unique<math::linal::GaussLinearSystemSolver>();
    }

    Vector coarse_solution(const Vector& rhs) const {
        return m_coarse_solver->solve(m_matrix, rhs);
    }

protected:
    Params m_params;
    
private:
    const AnyMatrix& m_matrix;
    std::unique_ptr<math::linal::ILinearSystemSolver> m_coarse_solver;
};

namespace {

    class BandImpl final : public math::linal::MGPreconditioner::Impl {
        using Matrix = math::linal::BandMatrix;

    public:
        BandImpl(const Matrix& matrix, const Params& params)
            : Impl(matrix, params)
        {
            m_levels.clear();
            m_levels.push_back(matrix);

            // Создаем иерархию уровней
            while (m_levels.back().get_width() > m_params.min_coarse_size && m_levels.size() < m_params.max_levels) {
                Matrix coarse = create_coarse_matrix(m_levels.back());
                m_levels.push_back(coarse);
            }

            coarse_sovler_init();
        }

    private:
        bool check_levels() const override {
            return !m_levels.empty() && m_levels.back().get_width() <= m_params.min_coarse_size;
        }

        void v_cycle(size_t level_idx, const Vector& rhs, Vector& x) const override {
            if (level_idx == m_levels.size() - 1) {
                x = coarse_solution(rhs);
                return;
            }

            const Matrix& A = m_levels[level_idx];

            // Предварительное сглаживание
            for (size_t i = 0; i < m_params.n_pre_smooth; ++i) {
                smooth(A, rhs, x, level_idx);
            }

            Vector residual = rhs - A * x;

            // Ограничение невязки на более грубый уровень
            Vector coarse_rhs = restrict(residual, level_idx);
            Vector coarse_correction(coarse_rhs.size(), 0.0);

            // Рекурсивный вызов для более грубого уровня
            v_cycle(level_idx + 1, coarse_rhs, coarse_correction);

            // Интерполяция коррекции на более тонкий уровень
            Vector correction = interpolate(coarse_correction, level_idx);
            x += correction;

            // Пост-сглаживание
            for (size_t i = 0; i < m_params.n_post_smooth; ++i) {
                smooth(A, rhs, x, level_idx);
            }
        }

        void w_cycle(size_t level_idx, const Vector& rhs, Vector& x) const override {
            if (level_idx == m_levels.size() - 1) {
                x = coarse_solution(rhs);
                return;
            }

            const Matrix& A = m_levels[level_idx];

            // Предварительное сглаживание
            for (size_t i = 0; i < m_params.n_pre_smooth; ++i) {
                smooth(A, rhs, x, level_idx);
            }

            Vector residual = rhs - A * x;

            // Ограничение невязки на более грубый уровень
            Vector coarse_rhs = restrict(residual, level_idx);
            Vector coarse_correction(coarse_rhs.size(), 0.0);

            // Два рекурсивных вызова для более грубого уровня (W-цикл)
            w_cycle(level_idx + 1, coarse_rhs, coarse_correction);
            w_cycle(level_idx + 1, coarse_rhs, coarse_correction);

            // Интерполяция коррекции на более тонкий уровень
            Vector correction = interpolate(coarse_correction, level_idx);
            x += correction;

            // Пост-сглаживание
            for (size_t i = 0; i < m_params.n_post_smooth; ++i) {
                smooth(A, rhs, x, level_idx);
            }
        }

        void f_cycle(size_t level_idx, const Vector& rhs, Vector& x) const override {
            if (level_idx == m_levels.size() - 1) {
                x = coarse_solution(rhs);
                return;
            }

            const Matrix& A = m_levels[level_idx];

            // Предварительное сглаживание
            for (size_t i = 0; i < m_params.n_pre_smooth; ++i) {
                smooth(A, rhs, x, level_idx);
            }

            Vector residual = rhs - A * x;

            // Ограничение невязки на более грубый уровень
            Vector coarse_rhs = restrict(residual, level_idx);
            Vector coarse_correction(coarse_rhs.size(), 0.0);

            // Рекурсивный вызов F-цикла для более грубого уровня
            f_cycle(level_idx + 1, coarse_rhs, coarse_correction);

            // Интерполяция коррекции на более тонкий уровень
            Vector correction = interpolate(coarse_correction, level_idx);
            x += correction;

            // Дополнительное сглаживание и V-цикл
            for (size_t i = 0; i < m_params.n_post_smooth; ++i) {
                smooth(A, rhs, x, level_idx);
            }

            residual = rhs - A * x;
            coarse_rhs = restrict(residual, level_idx);
            coarse_correction.assign(coarse_rhs.size(), 0.0);
            v_cycle(level_idx + 1, coarse_rhs, coarse_correction);
            correction = interpolate(coarse_correction, level_idx);
            x += correction;

            // Финальное пост-сглаживание
            for (size_t i = 0; i < m_params.n_post_smooth; ++i) {
                smooth(A, rhs, x, level_idx);
            }
        }

        Matrix create_coarse_matrix_algebraic(const Matrix& fine) {
            // Алгебраический многосеточный метод (AMG) - построение грубой матрицы
            size_t n = fine.get_width();

            // 1. Определение сильных связей и создание C/F разбиения
            std::vector<bool> is_c_point(n, false);
            std::vector<std::vector<size_t>> strong_connections(n);
            std::vector<size_t> lambda(n, 0); // Количество сильных связей

            // Вычисление сильных связей и лямбда (мера важности)
            for (size_t i = 0; i < n; ++i) {
                double max_strong = 0.0;
                double diag = std::abs(fine[i][0]);

                // Находим максимальный по модулю элемент в строке
                for (size_t j = 1; j < fine[i].size(); ++j) {
                    if (i + j < n) {
                        double val = std::abs(fine[i][j]);
                        if (val > max_strong) {
                            max_strong = val;
                        }
                    }
                }

                // Определяем сильные связи
                double threshold = m_params.strong_threshold * max_strong;
                for (size_t j = 1; j < fine[i].size(); ++j) {
                    if (i + j < n) {
                        double val = std::abs(fine[i][j]);
                        if (val >= threshold && val > 0) {
                            strong_connections[i].push_back(i + j);
                            lambda[i]++;
                        }
                    }
                }

                // Также учитываем симметричную часть (нижний треугольник)
                for (size_t j = 1; j < fine[i].size(); ++j) {
                    if (i >= j) {
                        double val = std::abs(fine[i - j][j]);
                        if (val >= threshold && val > 0) {
                            strong_connections[i].push_back(i - j);
                            lambda[i]++;
                        }
                    }
                }
            }

            // 2. Выбор C-точек (грубых точек) с использованием алгоритма Ruge-Stuben
            std::vector<size_t> c_points;
            std::vector<size_t> f_points;

            while (true) {
                // Находим точку с максимальным lambda
                size_t max_lambda = 0;
                size_t i_max = n;
                for (size_t i = 0; i < n; ++i) {
                    if (!is_c_point[i] && lambda[i] >= max_lambda) {
                        max_lambda = lambda[i];
                        i_max = i;
                    }
                }

                if (i_max == n) break; // Все точки обработаны

                // Добавляем как C-точку
                is_c_point[i_max] = true;
                c_points.push_back(i_max);

                // Обновляем lambda для соседей
                for (size_t j : strong_connections[i_max]) {
                    if (!is_c_point[j]) {
                        lambda[j]++; // Увеличиваем важность соседей

                        // Добавляем соседей сильных связей j как F-точки
                        for (size_t k : strong_connections[j]) {
                            if (!is_c_point[k] && k != i_max) {
                                is_c_point[k] = false;
                                lambda[k]--; // Уменьшаем важность
                            }
                        }
                    }
                }
            }

            // Остальные точки - F-точки
            for (size_t i = 0; i < n; ++i) {
                if (!is_c_point[i]) {
                    f_points.push_back(i);
                }
            }

            size_t coarse_n = c_points.size();
            if (coarse_n == 0) return Matrix();

            // 3. Построение интерполяционного оператора P и грубой матрицы (Galerkin projection)
            // A_coarse = P^T * A_fine * P

            // Для простоты используем прямую интерполяцию
            Matrix coarse(coarse_n);
            size_t coarse_isl = std::max<size_t>(1, fine.get_isl() / 2);
            coarse.set_isl(coarse_isl);

            // Создаем отображение индексов
            std::vector<size_t> index_map(n);
            for (size_t i = 0; i < coarse_n; ++i) {
                index_map[c_points[i]] = i;
            }

            // Построение грубой матрицы
            for (size_t i = 0; i < coarse_n; ++i) {
                size_t fine_i = c_points[i];

                for (size_t j = i; j < std::min(i + coarse_isl, coarse_n); ++j) {
                    size_t fine_j = c_points[j];

                    // Вычисляем элемент грубой матрицы
                    double sum = 0.0;

                    // Учитываем вклад от A_fine[i][k] * P[k][j]
                    // Для прямой интерполяции P - это просто отображение C-точек

                    // Диагональный элемент
                    if (fine_i == fine_j) {
                        sum += fine[fine_i][0];
                    }

                    // Внедиагональные элементы (верхний треугольник)
                    for (size_t k = 1; k < fine[fine_i].size(); ++k) {
                        size_t col = fine_i + k;
                        if (col < n && is_c_point[col]) {
                            sum += fine[fine_i][k] * (col == fine_j ? 1.0 : 0.0);
                        }
                    }

                    // Внедиагональные элементы (нижний треугольник - симметричная часть)
                    for (size_t k = 1; k < fine[fine_i].size(); ++k) {
                        size_t row = fine_i - k;
                        if (row < n && is_c_point[row]) {
                            sum += fine[row][k] * (row == fine_j ? 1.0 : 0.0);
                        }
                    }

                    coarse[i][j - i] = sum;
                }
            }

            return coarse;
        }

        Matrix create_coarse_matrix_geometric(const Matrix& fine) {
            // Здесь должна быть реализация создания грубой матрицы
            // В упрощенном виде (нынешнем) просто береться каждая вторая точкп

            size_t n = fine.get_width();
            size_t coarse_n = (n + 1) / 2;
            Matrix coarse(coarse_n);

            size_t coarse_isl = std::max<size_t>(1, fine.get_isl() / 2);
            coarse.set_isl(coarse_isl);

            for (size_t i = 0; i < coarse_n; ++i) {
                size_t fine_i = 2 * i;
                if (fine_i >= n) break;

                for (size_t j = i; j < std::min(i + coarse_isl, coarse_n); ++j) {
                    size_t fine_j = 2 * j;
                    if (fine_j >= n) break;

                    size_t band_offset = fine_j - fine_i;
                    if (band_offset < fine[fine_i].size()) {
                        coarse[i][j - i] = fine[fine_i][band_offset];
                    }
                }
            }

            return coarse;
        }

        Matrix create_coarse_matrix(const Matrix& fine) {
            switch (m_params.mg_type) {
                case MGType::Algebraic: return create_coarse_matrix_algebraic(fine);
                case MGType::Geometric: return create_coarse_matrix_geometric(fine);
            }
        }

        Vector restrict_algebraic(const Vector& fine, size_t level_idx) const override {
            const Matrix& fine_matrix = m_levels[level_idx];
            size_t n_fine = fine_matrix.get_width();
            size_t n_coarse = m_levels[level_idx + 1].get_width();

            // Если это первый уровень, используем инъекцию
            if (level_idx == 0) {
                Vector coarse(n_coarse, 0.0);
                for (size_t i = 0; i < n_coarse; ++i) {
                    size_t fine_i = 2 * i;
                    if (fine_i < n_fine) {
                        coarse[i] = fine[fine_i];
                    }
                }
                return coarse;
            }

            // Для более высоких уровней используем полное взвешенное ограничение
            Vector coarse(n_coarse, 0.0);
            const auto& fine_data = fine_matrix.data();

            for (size_t i = 0; i < n_coarse; ++i) {
                double sum_weights = 0.0;
                double sum_values = 0.0;

                // Учитываем диагональный элемент
                sum_weights += std::abs(fine_data[i][0]);
                sum_values += fine_data[i][0] * fine[i];

                // Учитываем соседей
                for (size_t j = 1; j < fine_data[i].size(); ++j) {
                    if (i + j < n_fine) {
                        double weight = std::abs(fine_data[i][j]);
                        sum_weights += weight;
                        sum_values += weight * fine[i + j];
                    }
                    if (i >= j) {
                        double weight = std::abs(fine_data[i - j][j]);
                        sum_weights += weight;
                        sum_values += weight * fine[i - j];
                    }
                }

                if (sum_weights > 0) {
                    coarse[i] = sum_values / sum_weights;
                }
            }

            return coarse;
        }

        Vector restrict_geometric(const Vector& fine, size_t level_idx) const override {
            // Простая инъекция - берем каждую вторую точку
            const Matrix& fine_matrix = m_levels[level_idx];
            size_t coarse_n = m_levels[level_idx + 1].get_width();
            Vector coarse(coarse_n, 0.0);

            for (size_t i = 0; i < coarse_n; ++i) {
                size_t fine_i = 2 * i;
                if (fine_i < fine_matrix.get_width()) {
                    coarse[i] = fine[fine_i];
                }
            }

            return coarse;
        }

        Vector interpolate(const Vector& coarse, size_t level_idx) const {
            // Линейная интерполяция
            size_t fine_size = m_levels[level_idx].get_width();
            Vector fine(fine_size, 0.0);

            for (size_t i = 0; i < coarse.size(); ++i) {
                if (2 * i < fine_size) {
                    fine[2 * i] = coarse[i];
                }
                if (2 * i + 1 < fine_size && i + 1 < coarse.size()) {
                    fine[2 * i + 1] = 0.5 * (coarse[i] + coarse[i + 1]);
                }
            }

            return fine;
        }

        void smooth(const Matrix& A, const Vector& rhs, Vector& x, size_t level_idx) const {
            switch (m_params.smoother_type) {
            case SmootherType::JACOBI:
                jacobi_smooth(A, rhs, x);
                break;
            case SmootherType::GAUSS_SEIDEL:
                gauss_seidel_smooth(A, rhs, x);
                break;
            case SmootherType::ILU0:
                ilu0_smooth(A, rhs, x, level_idx);
                break;
            }
        }

        void jacobi_smooth(const Matrix& A, const Vector& rhs, Vector& x) const {
            Vector new_x = x;
            size_t isl = A.get_isl();  
            size_t n = A.get_width();

            for (size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                double diag = A[i][0];  

                for (size_t j = 1; j <= isl && i + j < n; ++j) {
                    sum += A[i][j] * x[i + j];
                }

                for (size_t j = 1; j <= isl && i >= j; ++j) {
                    sum += A[i - j][j] * x[i - j];  
                }

                new_x[i] = (rhs[i] - sum) / diag;
            }

            x = new_x;
        }

        void gauss_seidel_smooth(const Matrix& A, const Vector& rhs, Vector& x) const {
            const size_t isl = A.get_isl();
            const size_t height = A.get_height();
            const size_t width = A.get_width();
            for (size_t i = 0; i < height; ++i) {
                double sum = 0.0;
                double diag = A[i][0]; 

                for (size_t j = 1; j < isl; ++j) {
                    size_t col = i + j;
                    if (col < width) {
                        sum += A[i][j] * x[col];
                    }
                }

                for (size_t j = 1; j < isl; ++j) {
                    size_t row = i - j;
                    if (row < height) {
                        sum += A[row][j] * x[row];
                    }
                }

                x[i] = (rhs[i] - sum) / diag;
            }
        }

        void ilu0_smooth(const Matrix& A, const Vector& rhs, Vector& x, size_t level_idx) const {
            // 1. Разложение ILU(0)
            size_t n = A.get_width();
            size_t isl = A.get_isl();

            // Создаем копии L и U частей
            std::vector<Vector> L(n), U(n);

            for (size_t i = 0; i < n; ++i) {
                L[i].resize(isl + 1, 0.0);
                U[i].resize(isl + 1, 0.0);

                // Инициализация L и U
                for (size_t j = 0; j <= isl && i + j < n; ++j) {
                    U[i][j] = A[i][j];
                }
                if (i < n) {
                    L[i][0] = 1.0; // Диагональ L
                }
            }

            // Выполняем неполное разложение
            for (size_t i = 0; i < n; ++i) {
                // Обрабатываем только ненулевые элементы в пределах полосы
                for (size_t k = 0; k < i && k <= isl; ++k) {
                    if (U[i - k - 1][k + 1] != 0.0 && U[i][0] != 0.0) {
                        double factor = U[i - k - 1][k + 1] / U[i][0];
                        L[i][k + 1] = factor;

                        // Обновляем строку i
                        for (size_t j = 0; j <= isl && i + j < n; ++j) {
                            if (U[i][j] != 0.0) {
                                size_t row = i - k - 1;
                                size_t col = j + k + 1;
                                if (row < n && col < n && col - row <= isl) {
                                    U[row][col - row] -= factor * U[i][j];
                                }
                            }
                        }
                    }
                }
            }

            // Ly = rhs 
            Vector y(n, 0.0);
            for (size_t i = 0; i < n; ++i) {
                double sum = rhs[i];
                for (size_t j = 1; j <= isl && i >= j; ++j) {
                    sum -= L[i][j] * y[i - j];
                }
                y[i] = sum / L[i][0];
            }

            // Ux = y 
            for (int i = n - 1; i >= 0; --i) {
                double sum = y[i];
                for (size_t j = 1; j <= isl && i + j < n; ++j) {
                    sum -= U[i][j] * x[i + j];
                }
                x[i] = sum / U[i][0];
            }
        }

    private:
        Params m_params;
        std::vector<Matrix> m_levels;
    };

    class DenseImpl final : public math::linal::MGPreconditioner::Impl {
        using Matrix = math::linal::DenseMatrix;

    public:
        DenseImpl(const Matrix& matrix, const Params& params)
            : Impl(matrix, params)
        {
            // ...
        }

    private:
        bool check_levels() const override {
            return true;
        }

        void v_cycle(size_t level_idx, const Vector& rhs, Vector& x) const override {
            // ...
        }

        void w_cycle(size_t level_idx, const Vector& rhs, Vector& x) const override {
            // ...
        }

        void f_cycle(size_t level_idx, const Vector& rhs, Vector& x) const override {
            // ...
        }

        Vector restrict_algebraic(const Vector& fine, size_t level_idx) const override {
            // ...
            return {};
        }

        Vector restrict_geometric(const Vector& fine, size_t level_idx) const override {
            // ...
            return {};
        }
    };

    class SparseImpl final : public math::linal::MGPreconditioner::Impl {
        using Matrix = math::linal::SparseMatrix;

    public:
        SparseImpl(const Matrix& matrix, const Params& params)
            : Impl(matrix, params)
        {
            // ...
        }

    private:
        bool check_levels() const override {
            return true;
        }

        void v_cycle(size_t level_idx, const Vector& rhs, Vector& x) const override {
            // ...
        }

        void w_cycle(size_t level_idx, const Vector& rhs, Vector& x) const override {
            // ...
        }

        void f_cycle(size_t level_idx, const Vector& rhs, Vector& x) const override {
            // ...
        }

        Vector restrict_algebraic(const Vector& fine, size_t level_idx) const override {
            // ...
            return {};
        }

        Vector restrict_geometric(const Vector& fine, size_t level_idx) const override {
            // ...
            return {};
        }
    };

} // namespace

math::linal::MGPreconditioner::MGPreconditioner(const Params& params)
    : m_params(params) {
}

void math::linal::MGPreconditioner::init(const AnyMatrix& matrix) {
    if (std::holds_alternative<BandMatrix>(matrix)) {
        m_impl = std::make_unique<BandImpl>(std::get<BandMatrix>(matrix), m_params);
    }
    else if (std::holds_alternative<DenseMatrix>(matrix)) {
        m_impl = std::make_unique<DenseImpl>(std::get<DenseMatrix>(matrix), m_params);
    }
    else if (std::holds_alternative<SparseMatrix>(matrix)) {
        m_impl = std::make_unique<SparseImpl>(std::get<SparseMatrix>(matrix), m_params);
    }
    else {
        throw std::runtime_error("Unsupported matrix type for MG preconditioner");
    }
}
