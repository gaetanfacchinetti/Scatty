#ifndef INTEGRATION_H
#define INTEGRATION_H

#define EPSILON_ROOTS 1e-14  // Precision for Newton's method


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>


void legendre(int n, double x, double* P, double* dP);
int gauss_legendre_points_and_weights(int n, double* x, double* w);

#endif // INTEGRATION_H