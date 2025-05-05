#include "GMRESSystemSolver.h"

#include <stdexcept>

math::linal::GMRESLinearSystemSolver::GMRESLinearSystemSolver(size_t Krylov_subspace_dimension, const IterativeSolvingParams& iter_params, const Params& params, std::unique_ptr<IPreconditioner> preconditioner)
    : IKrylovTypeLinearSystemSolver(Krylov_subspace_dimension, iter_params, params, std::move(preconditioner))
{
}

math::linal::FVector math::linal::GMRESLinearSystemSolver::solve(const AnyMatrix& matrix, const FVector& rhs, const FVector& x0) {
    return std::visit([this, &rhs, &x0](const auto& matrix) -> FVector {
        auto [need_to_solve, x, rhs_norm, r, beta, resid] = init_method(matrix, rhs, x0);
        if (!need_to_solve) {
            return x;
        }

        const size_t n = rhs.size();
        const size_t m = get_Krylov_subspace_dimension(n);
        const size_t max_iteration_count = get_max_iteration_count(n);

        std::vector<FVector> V(m + 1, FVector(n));  // ����� �������
        std::vector<double> H((m + 1) * m, 0.0);      // ������� �����������
        std::vector<double> cs(m + 1);              // �������� ��������
        std::vector<double> sn(m + 1);              // ������ ��������
        std::vector<double> s(m + 1);               // ������ ������ �����
        FVector xh(n);                            // ��������� ������

        for (size_t j = 1; j <= max_iteration_count;) {
            V[0] = r / beta; // ������ �������� ������
            s[0] = beta;
            std::fill(s.begin() + 1, s.end(), 0.0);

            for (size_t i = 0; i < m && j <= max_iteration_count; ++i, ++j) {
                xh = apply_precondition(V[i]);

                FVector w = matrix * xh;

                // ��������������� (���������������� �����-�����)
                for (size_t k = 0; k <= i; ++k) {
                    H[k + i * (m + 1)] = w.dot(V[k]);
                    w -= V[k] * H[k + i * (m + 1)];
                }

                // ������������
                H[(i + 1) + i * (m + 1)] = w.norm();
                if (H[(i + 1) + i * (m + 1)] < m_iterative_solving_params.tolerance) {
                    // ��������� ������� ��������� ������
                    break;
                }
                if (i < m) {
                    V[i + 1] = w / H[(i + 1) + i * (m + 1)];
                }

                // ���������� ���������� �������� � ������ ������� H
                for (size_t k = 0; k < i; ++k) {
                    apply_Givens_rotation(H[k + i * (m + 1)], H[k + 1 + i * (m + 1)], cs[k], sn[k]);
                }

                std::tie(cs[i], sn[i]) = generate_Givens_rotation(H[i + i * (m + 1)], H[i + 1 + i * (m + 1)]);

                // ���������� �������� � H � s
                apply_Givens_rotation(H[i + i * (m + 1)], H[i + 1 + i * (m + 1)], cs[i], sn[i]);
                apply_Givens_rotation(s[i], s[i + 1], cs[i], sn[i]);

                resid = fabs(s[i + 1] / rhs_norm);
                if (resid < m_iterative_solving_params.tolerance) {
                    compute_correction(i + 1, n, H, m + 1, s, V, x);
                    collect_solution_stats({ resid, j });
                    return x;
                }
            }

            // ���������� ������� ����� m ��������
            compute_correction(m, n, H, m + 1, s, V, x);

            r = rhs - matrix * x;
            beta = r.norm();
            resid = beta / rhs_norm;

            if (check_convergence_criterion(r, rhs_norm, j)) {
                collect_solution_stats({ resid, j });
                return x;
            }
        }

        collect_solution_stats({ resid, max_iteration_count });
        return x;
        }, matrix);
}

// ������� ������� ����������� ������� H y = s � ���������� x += Z y
void math::linal::GMRESLinearSystemSolver::compute_correction(size_t k, size_t n, const std::vector<double>& H, size_t ldH, const std::vector<double>& s, const std::vector<FVector>& V, FVector& x) const {
    auto y = solve_H_system(k, H, ldH, s, m_params.throw_exceptions);

    // ���������� ���������: x += M*V*y
    FVector correction(n, 0.0);
    for (size_t i = 0; i < k; ++i) {
        correction += V[i] * y[i];
    }

    x += apply_precondition(correction);
}

math::linal::FVector math::linal::GMRESLinearSystemSolver::solve_H_system(size_t k, const std::vector<double>& H, size_t ldH, const std::vector<double>& s, bool throw_exception) {
    FVector y(k);
    for (size_t i = 0; i < k; ++i) {
        y[i] = s[i];
    }

    // �������� ��� ��� ������� ����������� �������
    for (int i = k - 1; i >= 0; --i) {
        double diag = H[i + i * ldH];

        // ������ �� ������� �� ����
        if (dcmp(diag) == 0) {
            if (dcmp(y[i]) == 0) {
                y[i] = 0.0;
            }
            else {
                if (throw_exception)
                    throw std::runtime_error("GMRES: singular H matrix encountered");
                break;
            }
        }
        else {
            y[i] /= diag;
        }

        for (int j = 0; j < i; ++j) {
            y[j] -= H[j + i * ldH] * y[i];
        }
    }

    return y;
}

// ��������������� ������� ��� �������� �������
std::pair<double, double> math::linal::GMRESLinearSystemSolver::generate_Givens_rotation(double dx, double dy) {
    if (math::dcmp(dy) == 0) {
        return { 1.0, 0.0 }; // cs, sn
    }
    else if (fabs(dy) > fabs(dx)) {
        double tmp = dx / dy;
        double sn = 1.0 / sqrt(1.0 + tmp * tmp);
        return { tmp * sn, sn };
    }
    else {
        double tmp = dy / dx;
        double cs = 1.0 / sqrt(1.0 + tmp * tmp);
        return { cs, tmp * cs };
    }
}

void math::linal::GMRESLinearSystemSolver::apply_Givens_rotation(double& dx, double& dy, double cs, double sn) {
    if (std::isnan(cs) || std::isnan(sn)) {
        dx = dy = 0.0;
        return;
    }
    const double tmp = cs * dx + sn * dy;
    dy = cs * dy - sn * dx;
    dx = tmp;
}