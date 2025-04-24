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

#include "integration.h"


// Function to compute the Legendre polynomial P_n(x) and its derivative
void legendre(int32_t n, double x, double* P, double* dP) {
    
    double P0 = 1.0, P1 = x, Pk = 0;
    double dP1 = 1.0, dPk = 0;
    
    for (int k = 2; k <= n; k++) {

        Pk = ((2.0 * k - 1.0) * x * P1 - (k - 1.0) * P0) / k;
        dPk = k * (P1 - x * Pk) / (1.0 - x * x);
        
        P0 = P1;
        P1 = Pk;
        // dP0 = dP1;
        dP1 = dPk;
    }
    
    *P = P1;
    *dP = dP1;
}


// Function to compute Gauss-Legendre quadrature points and weights
int set_gauss_legendre_points_and_weights(int32_t n, double* x, double* w) {
    
    int m = (n + 1) / 2; // Only need to compute half the roots due to symmetry
    double xi = 0, P = 0, dP = 0, dx=0;
    
    // define a maximal bound for the newton solver
    int max_iter = 1000000, iter=0;

    for (int32_t i = 0; i < m; i++) {
        // Initial guess for the root using an approximation
        xi = cos(M_PI * (i + 0.75) / (n + 0.5));

        // Newton's method to find the root
        do {
            legendre(n, xi, &P, &dP);
            dx = -P / dP;
            xi += dx;
            iter = iter+1;
            
        } while ((fabs(dx) > EPSILON_ROOTS) && iter < max_iter);
        
        if (iter == max_iter) return -1;

        // save the computed root
        x[i] = -xi;
        x[n - i - 1] = xi;

        // compute the corresponding weight
        w[i] = w[n - i - 1] = 2.0 / ((1.0 - xi * xi) * dP * dP);
    }
    
    return 0;
}


// Function to compute Gauss-Legendre quadrature points and weights
int set_root_and_weights_scale(int32_t n, double a, double b, double* x, double* w, bool in_log) {
    
    if (in_log == true){
        
        // first check that the boundaries make sense
        if ((a < 0) || (b < 0)) {
            return -1;
        }
        
        // redefine the roots and weights
        for (int32_t i = 0; i < n; i++){
            x[i] = pow(b, (1.0 + x[i])/2.0) * pow(a, (1.0 - x[i])/2.0);
            w[i] = x[i] * log(b/a)/2.0 * w[i];
            
        }
    } else {

        // redefine the roots and weights
        for (int32_t i = 0; i < n; i++){
            w[i] = (b-a)/2.0 * w[i];
            x[i] = ((b-a) * x[i] + a + b)/2.0;
        }
    }

    return 0;
}


// uses the Box-Muller transformation
void draw_gauss_monte_carlo_points(int32_t n, double* x, double mu, double sigma)
{
    // set the seed for the rand() function
    srand(time(0));

    double u, v;

    for (int32_t i = 0; i < n; i = i+2){
        
        u = (double) rand() / RAND_MAX;
        v = (double) rand() / RAND_MAX;
    
        x[i] = mu + sigma * sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
        
        if (i+1 < n){
            x[i+1] = mu + sigma * sqrt(-2.0 * log(u)) * sin(2.0 * M_PI * v);
        }
    }

}


// Test functions for the implementation of the integral
// Follow this example to compute your own integral

double test_function(double x){
    return 1/x + 1/(x*x);
}

double test_integral(){
    
    int err_0, err_1, n = 10;

    // set the weights and roots
    double x[n], w[n];
    err_0 = set_gauss_legendre_points_and_weights(n, x, w);
    err_1 = set_root_and_weights_scale(n, 0.01, 100.0, x, w, true);

    if (err_0 < 0 || err_1 < 0){
        return -99.0;
    }

    double res = 0;
    for (int32_t i = 0; i < n; i++){
        res = res + test_function(x[i]) * w[i];
    }

    return res;

}