#include <Math/Calculus/Integration/Integrator.h>

double math::calcls::Integrator::trapezoidal(const std::function<double(double)>& f, double a, double b, int n) {
    if (n <= 0) return 0.0;

    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(x);
    }

    return sum * h;
}

double math::calcls::Integrator::simpson(const std::function<double(double)>& f, double a, double b, int n) {
    if (n <= 0 || n % 2 != 0) return 0.0; // n должно быть четным

    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i += 2) {
        double x = a + i * h;
        sum += 4.0 * f(x);
    }

    for (int i = 2; i < n; i += 2) {
        double x = a + i * h;
        sum += 2.0 * f(x);
    }

    return sum * h / 3.0;
}
