/** ----------------------------------------------------------------------
 * This file is part of Scatty.
 *
 * Copyright (c) 2024, Gaétan Facchinetti
 *
 * Scatty is free software: you can redistribute it and/or modify it 
 * under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or any 
 * later version. Scatty is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU 
 * General Public License along with NNERO. 
 * If not, see <https://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------
 */

#ifndef INTEGRATION_H
#define INTEGRATION_H

// Precision for Newton's method
#define EPSILON_ROOTS 1e-14  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>


/**
 * @brief Legendre polynomial of degree n and its derivative evaluated at x
 *
 * @param n    Degree of polynomial (int32_t*).
 * @param x    Point of evaluation (double).
 * @param P    Value of the polynomial (double*).
 * @param dP   Value of the polynomial derivative (double*).
 * @return     nothing
 */
void legendre(int32_t n, double x, double* P, double* dP);

/**
 * @brief Sets nodes and weights of Gauss-Legendre integration routine
 *        in the domain [-1, 1]
 *
 * @param n    Number of nodes (int32_t*).
 * @param x    Values of nodes (double*).
 * @param P    Weights of nodes (double*).
 * @return     0 if the routine exited normally (int).
 */
int set_gauss_legendre_points_and_weights(int32_t n, double* x, double* w);

/**
 * @brief Rescales the nodes and weights for an interval of integration [a, b]
 *        and with log scaling or not
 *
 * @param n       Number of nodes (int32_t*).
 * @param a       Min bound(double).
 * @param b       Max bound (double).
 * @param x       Values of nodes for the domain [-1, 1] (double*).
 * @param w       Values of weights for the domain [-1, 1] (double*).
 * @param in_log  Set to true for log scaling (bool).
 * @return        0 if the routine exited normally (int).
 */
int set_root_and_weights_scale(int32_t n, double a, double b, double* x, double* w, bool in_log);

/**
 * @brief Draw n points from a Gaussian distribution
 *
 * @param n       Number of nodes (int32_t*).
 * @param x       Values of nodes (double*).
 * @param mu      Average of the Gaussian (double).
 * @param sigma   Standard deviation of the Gaussian (double).
 * @return        0 if the routine exited normally (int).
 */
void draw_gauss_monte_carlo_points(int32_t n, double* x, double mu, double sigma);


// Test functions
double test_function(double x);
double test_integral();

#endif // INTEGRATION_H