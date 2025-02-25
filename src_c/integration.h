#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <gsl/gsl_complex.h>

double integrate_trapezoidal(double (*func)(double), double a, double b, int n);


#endif // INTEGRATION_H