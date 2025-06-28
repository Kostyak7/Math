#include <Math/Calculus/ODESolvers/EulerSolver.h>

#include <Math/Common.h>

double math::calcls::EulerSolver::solve(const std::function<double(double, double)>& f, double y0, double x0, double x, double h = 1e-6) override {
    if (h <= 0.0) return y0;

    double y = y0;

    for (double current_x = x0; std::abs(current_x - x) > TOLERANCE; current_x += h) {
        if (current_x + h > x) {
            h = x - current_x;
        }

        y += h * f(current_x, y);
    }

    return y;
}