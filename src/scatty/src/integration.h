#ifndef INTEGRATION_H
#define INTEGRATION_H

#define EPSILON_ROOTS 1e-14  // Precision for Newton's method


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>


void legendre(int n, double x, double* P, double* dP);
int gauss_legendre_points_and_weights(int n, double* x, double* w);
int set_root_and_weights_scale(int n, double a, double b, double* x, double* w, bool in_log);

double test_function();
double test_integral();

#endif // INTEGRATION_H