#include "integration.h"

double integrate_trapezoidal(double (*func)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (func(a) + func(b));

    for (int i = 1; i < n; i++) {
        sum += func(a + i * h);
    }

    return sum * h;
}



