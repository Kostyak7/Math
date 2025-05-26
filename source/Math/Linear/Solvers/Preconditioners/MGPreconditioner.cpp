#include <Utils/Math/Linear/Solvers/Preconditioners/MGPreconditioner.h>

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

    using Vector = DVector;

public:
    Impl(const Params& params)
        : m_params(params)
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

    Vector coarse_solution(const AnyMatrixConstRef& matrix, const Vector& rhs) const {
        return m_coarse_solver->solve(matrix, rhs);
    }

protected:
    Params m_params;
    
private:
    std::unique_ptr<math::linal::ILinearSystemSolver> m_coarse_solver;
};

namespace {

    class BandImpl final : public math::linal::MGPreconditioner::Impl {
        using Matrix = math::linal::BandMatrix;

    public:
        BandImpl(const Matrix& matrix, const Params& params)
            : Impl(params)
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
                x = coarse_solution(m_levels.front(), rhs);
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
                x = coarse_solution(m_levels.front(), rhs);
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
                x = coarse_solution(m_levels.front(), rhs);
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
            size_t n = fine.get_width();

            // 1. Улучшенный алгоритм выбора C-точек с учетом симметрии
            std::vector<bool> is_c_point(n, false);
            std::vector<double> weights(n, 0.0);

            // Вычисляем веса для точек
            for (size_t i = 0; i < n; ++i) {
                double diag = std::abs(fine[i][0]);
                double sum = 0.0;
                for (size_t j = 1; j < fine[i].size(); ++j) {
                    if (i + j < n) sum += std::abs(fine[i][j]);
                    if (i >= j) sum += std::abs(fine[i - j][j]);
                }
                weights[i] = diag / (sum + 1e-12);
            }

            // Выбираем C-точки
            std::vector<size_t> c_points;
            for (size_t i = 0; i < n; ++i) {
                if (weights[i] > m_params.strong_threshold) {
                    is_c_point[i] = true;
                    c_points.push_back(i);
                }
            }

            size_t coarse_n = c_points.size();
            if (coarse_n == 0) return Matrix();

            // 2. Построение грубой матрицы с учетом симметрии
            Matrix coarse(coarse_n);
            size_t coarse_isl = std::max<size_t>(1, fine.get_isl() / 2);
            coarse.set_isl(coarse_isl);

            for (size_t i = 0; i < coarse_n; ++i) {
                size_t fine_i = c_points[i];
                for (size_t j = i; j < std::min(i + coarse_isl, coarse_n); ++j) {
                    size_t fine_j = c_points[j];

                    // Учитываем симметричные элементы
                    double val = 0.0;
                    if (fine_i == fine_j) {
                        val = fine[fine_i][0];
                    }
                    else if (fine_j > fine_i && (fine_j - fine_i) < fine[fine_i].size()) {
                        val = fine[fine_i][fine_j - fine_i];
                    }
                    else if (fine_i > fine_j && (fine_i - fine_j) < fine[fine_j].size()) {
                        val = fine[fine_j][fine_i - fine_j];
                    }

                    coarse[i][j - i] = val;
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

            Vector coarse(n_coarse, 0.0);

            // Для AMG лучше использовать полное взвешенное ограничение на всех уровнях
            for (size_t i = 0; i < n_coarse; ++i) {
                double sum_weights = 0.0;
                double sum_values = 0.0;

                // Учитываем диагональный элемент
                sum_weights += std::abs(fine_matrix[i][0]);
                sum_values += fine_matrix[i][0] * fine[i];

                // Учитываем соседей (верхний треугольник)
                for (size_t j = 1; j < fine_matrix[i].size(); ++j) {
                    if (i + j < n_fine) {
                        double weight = std::abs(fine_matrix[i][j]);
                        sum_weights += weight;
                        sum_values += weight * fine[i + j];
                    }
                }

                // Учитываем соседей (нижний треугольник - симметричная часть)
                for (size_t j = 1; j < fine_matrix[i].size(); ++j) {
                    if (i >= j) {
                        double weight = std::abs(fine_matrix[i - j][j]);
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

            // Для AMG лучше использовать более сложную интерполяцию
            if (level_idx > 0) {
                // Используем информацию о связях из матрицы
                const Matrix& A = m_levels[level_idx - 1];
                for (size_t i = 0; i < fine_size; ++i) {
                    double sum = 0.0;
                    double sum_weights = 0.0;

                    // Учитываем диагональный элемент
                    sum += A[i][0] * coarse[i];
                    sum_weights += std::abs(A[i][0]);

                    // Учитываем соседей
                    for (size_t j = 1; j < A[i].size(); ++j) {
                        if (i + j < fine_size) {
                            sum += A[i][j] * coarse[i + j];
                            sum_weights += std::abs(A[i][j]);
                        }
                        if (i >= j) {
                            sum += A[i - j][j] * coarse[i - j];
                            sum_weights += std::abs(A[i - j][j]);
                        }
                    }

                    if (sum_weights > 0) {
                        fine[i] = sum / sum_weights;
                    }
                }
            }
            else {
                // Для самого грубого уровня - линейная интерполяция
                for (size_t i = 0; i < coarse.size(); ++i) {
                    if (2 * i < fine_size) fine[2 * i] = coarse[i];
                    if (2 * i + 1 < fine_size) {
                        fine[2 * i + 1] = 0.5 * (coarse[i] + (i + 1 < coarse.size() ? coarse[i + 1] : 0));
                    }
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
            size_t n = A.get_width();
            size_t isl = A.get_isl();

            // Временные структуры для L и U
            std::vector<std::unordered_map<size_t, double>> L(n), U(n);

            // Инициализация
            for (size_t i = 0; i < n; ++i) {
                // Диагональ U
                U[i][i] = A[i][0];

                // Внедиагональные элементы U
                for (size_t j = 1; j < A[i].size() && i + j < n; ++j) {
                    U[i][i + j] = A[i][j];
                }

                // Диагональ L
                L[i][i] = 1.0;
            }

            // Неполное разложение
            for (size_t i = 0; i < n; ++i) {
                // Обработка только в пределах полосы
                for (size_t k = 0; k < i && k <= isl; ++k) {
                    size_t j = i - k - 1;
                    if (U[j].count(i) && U[i].count(i)) {
                        double factor = U[j][i] / U[i][i];
                        L[i][j] = factor;

                        // Обновление строки j
                        for (const auto& elem : U[i]) {
                            size_t col = elem.first;
                            if (col >= i && col - j <= isl) { // В пределах полосы
                                U[j][col] -= factor * elem.second;
                            }
                        }
                    }
                }
            }

            // Решение Ly = rhs
            Vector y(n, 0.0);
            for (size_t i = 0; i < n; ++i) {
                double sum = rhs[i];
                for (const auto& elem : L[i]) {
                    size_t j = elem.first;
                    if (j < i) sum -= elem.second * y[j];
                }
                y[i] = sum / L[i][i];
            }

            // Решение Ux = y
            for (int i = n - 1; i >= 0; --i) {
                double sum = y[i];
                for (const auto& elem : U[i]) {
                    size_t j = elem.first;
                    if (j > i) sum -= elem.second * x[j];
                }
                x[i] = sum / U[i][i];
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
            : Impl(params)
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
            : Impl(params)
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

void math::linal::MGPreconditioner::init(const AnyMatrixConstRef& matrix) {
    if (std::holds_alternative<std::reference_wrapper<const BandMatrix>>(matrix)) {
        m_impl = std::make_unique<BandImpl>(std::get<std::reference_wrapper<const BandMatrix>>(matrix).get(), m_params);
    }
    else if (std::holds_alternative<std::reference_wrapper<const DenseMatrix>>(matrix)) {
        m_impl = std::make_unique<DenseImpl>(std::get<std::reference_wrapper<const DenseMatrix>>(matrix).get(), m_params);
    }
    else if (std::holds_alternative<std::reference_wrapper<const SparseMatrix>>(matrix)) {
        m_impl = std::make_unique<SparseImpl>(std::get<std::reference_wrapper<const SparseMatrix>>(matrix).get(), m_params);
    }
    else {
        throw std::runtime_error("Unsupported matrix type for MG preconditioner");
    }
}
