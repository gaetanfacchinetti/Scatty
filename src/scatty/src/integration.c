#include "integration.h"


// Function to compute the Legendre polynomial P_n(x) and its derivative
void legendre(int n, double x, double* P, double* dP) {
    
    double P0 = 1.0, P1 = x, Pk = 0;
    double dP0 = 0.0, dP1 = 1.0, dPk = 0;
    
    for (int k = 2; k <= n; k++) {

        Pk = ((2.0 * k - 1.0) * x * P1 - (k - 1.0) * P0) / k;
        dPk = k * (P1 - x * Pk) / (1.0 - x * x);
        
        P0 = P1;
        P1 = Pk;
        dP0 = dP1;
        dP1 = dPk;
    }
    
    *P = P1;
    *dP = dP1;
}


// Function to compute Gauss-Legendre quadrature points and weights
int set_gauss_legendre_points_and_weights(int n, double* x, double* w) {
    
    int m = (n + 1) / 2; // Only need to compute half the roots due to symmetry
    double xi = 0, P = 0, dP = 0, dx=0;
    
    // define a maximal bound for the newton solver
    int max_iter = 1000000, iter=0;


    for (int i = 0; i < m; i++) {
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
int set_root_and_weights_scale(int n, double a, double b, double* x, double* w, bool log) {
    
    int m = (n + 1) / 2;
    
    double yi = 0;

    //double 

    for (int i = 0; i < m; i++){
        
        yi = ((b-a)* x[i] + a + b)/2.0;
        x[i] = yi;
        x[n - i - 1] = -yi;

        w[i] = w[n - i - 1] = (b-a)/2.0 * w[i];
    }

    return 0;
}
