#include <Math/Calculus/Differentiation.h>

double MATH_EXPORT derivative(const std::function<double(double)>& f, double x, double h = 1e-5) {
    return (f(x + h) - f(x - h)) / (2 * h);
}