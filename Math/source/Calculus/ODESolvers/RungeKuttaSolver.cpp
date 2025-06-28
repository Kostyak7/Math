#include <Math/Calculus/ODESolvers/RungeKuttaSolver.h>

#include <Math/Common.h>

double math::calcls::RungeKuttaSolver::solve(const std::function<double(double, double)>& f, double y0, double x0, double x, double h = 1e-6) override {
    if (h <= 0.0) return y0;

    double y = y0;

    for (double current_x = x0; std::abs(current_x - x) > TOLERANCE; current_x += h) {
        if (current_x + h > x) {
            h = x - current_x;
        }

        double k1 = h * f(current_x, y);
        double k2 = h * f(current_x + h / 2, y + k1 / 2);
        double k3 = h * f(current_x + h / 2, y + k2 / 2);
        double k4 = h * f(current_x + h, y + k3);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }

    return y;
}