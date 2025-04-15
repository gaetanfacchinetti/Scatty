#include "cross_section.h"

double coulomb_phase_shift(int l, double zeta)
{

    // convention
    if (l == -1) return 0.0;

    // first compute delta for l = 0
    double delta_0 = - zeta * GAMMA_E;
    
    double new_term = 0;
    int m = 0;

    do {
        new_term = ( atan(zeta/(m+1.0)) - zeta / (m+1.0) );
        delta_0 = delta_0 - new_term;
        m = m+1;

    } while (m < 10000 && fabs(new_term/delta_0) > 1e-3);
    
    if (l == 0)
        return delta_0;

    // if l > 0 need to add more terms
    double delta_l = delta_0;

    for (m=1; m <=l; m++)
        delta_l = delta_l + atan(zeta/(1.0*m));

    return delta_l;
}


double* coulomb_phase_shift_arr(int l, double* zeta, int zeta_size)
{
    double* delta_l = (double*)malloc(zeta_size * sizeof(double));
    if (!delta_l) return NULL; // handle memory allocation failure

    // first get the values for the minimum l
    // assumes the values of l are given in ascending order
    for (int j = 0; j < zeta_size ; j++) {
        delta_l[j] = coulomb_phase_shift(l, zeta[j]);
    }
            
    return delta_l;
}


double* coulomb_phase_shift_grid(int* l, double* zeta, int l_size, int zeta_size)
{
    double* delta_l = (double*)malloc(l_size * zeta_size * sizeof(double));
    if (!delta_l) return NULL; // handle memory allocation failure

    int i_init = 0;

    // special case where the series starts with lm1
    if (l[0] == -1) {
        i_init = 1;
        for (int j = 0; j < zeta_size; j++) {
            delta_l[j] = 0.0;
        }
    }
        
    // first get the values for the minimum l
    // assumes the values of l are given in ascending order
    for (int j = 0; j < zeta_size; j++) {
        delta_l[i_init * zeta_size + j] = coulomb_phase_shift(l[i_init], zeta[j]);
    }
  
    // then keeps increasing for highest values of l
    // values of l must be given as a contiguous list
    for (int i = i_init + 1; i < l_size; i++) {
        for (int j = 0; j < zeta_size; j++){
            delta_l[i * zeta_size + j] = delta_l[(i-1)*zeta_size + j] +  atan(zeta[j]/(1.0*l[i]));
        }
    }

    return delta_l;
}


double transfer_factor(int l, double delta_l, double delta_lm1) {
    return sin(delta_l) * ( (1+l) * sin(delta_l) + l * sin(delta_l - 2*delta_lm1) );
}

double n_coulomb_transfer_cross_section(int l, double zeta)
{
    // only compute the phase shift once and get the adjacent values fast
    double delta_lm1 = l > 0 ? coulomb_phase_shift(l-1, zeta) : 0.0;
    double delta_l   = l > 0 ? delta_lm1 + atan(zeta/(1.0*l)) : coulomb_phase_shift(0, zeta);
    
    return transfer_factor(l, delta_l, delta_lm1) / zeta / zeta;
}


double* n_coulomb_transfer_cross_section_grid(int* l, double* zeta, int l_size, int zeta_size)
{
    // first get the values of the phase shifts
    int* lm1 = (int*)malloc((l_size+1) * sizeof(int));
    if (!lm1) return NULL;

    lm1[0] = l[0]-1;
    memmove(lm1 + 1, l, l_size * sizeof(int));

    double* delta_lm1 = coulomb_phase_shift_grid(lm1, zeta, l_size+1, zeta_size);
    
    // then compute the array of normalised l-wave cross-section 
    double* s_l = (double*)malloc(l_size * zeta_size * sizeof(double));
    if (!s_l) return NULL;

    for (int i = 0; i < l_size; i++) {
        for (int j = 0; j < zeta_size; j++){
            s_l[i * zeta_size + j] = transfer_factor(l[i], delta_lm1[(i+1) * zeta_size + j], delta_lm1[i * zeta_size + j]) / zeta[j] / zeta[j];
        }
    }

    // free the arrays of ls and deltas
    free_int_ptr(lm1);
    free_double_ptr(delta_lm1);

    return s_l;
}


double* n_coulomb_ur_transfer_cross_section_arr(int*l, int l_size)
{
    double* s_l = (double*)malloc(l_size * sizeof(double));
    if (!s_l) return NULL;

    s_l[0] = pow(GAMMA_E, 2);

    // initalise the digamma function at 1 (l=0, l+1=1)
    double psi = - GAMMA_E;

    for (int i = 1; i < l_size; i++){
        // psi(1+l) = psi(l) + 1/l
        psi = psi + 1/l[i];
        s_l[i] = psi*(psi+2);
    } 

    return s_l;
}



double mu_I(double y)
{
    return ((y-1) + exp(-2*y)*(y+1)) / (y*y);
}



double* r_chi_coulomb_arr(int* l, double zeta_r, double w_r, int l_size, int n)
{
    // --------------------------
    // set the weights and roots
    int err_0, err_1;
    double x[n], w[n];
    err_0 = set_gauss_legendre_points_and_weights(n, x, w);
    err_1 = set_root_and_weights_scale(n, fmax(w_r*(-5.0 + w_r), 0.0), w_r*(5.0 + w_r), x, w, false);
    if (err_0 < 0 || err_1 < 0) return NULL;


    // --------------------------
    // prepare the integration

    // get the grid of values on which to evaluate the cross-section
    double* zeta = (double*)malloc(n*sizeof(double));
    if (!zeta) return NULL;

    for (int j=0; j< n; j++) 
        zeta[j] = w_r * zeta_r / x[j];

    // get the values of s_l we need to integrate over
    double* s_l = n_coulomb_transfer_cross_section_grid(l, zeta, l_size, n);

    // preparing the returned array
    double* r_l = (double*)malloc(l_size * sizeof(double));
    if (!r_l) return NULL;

    // ----------------------------
    // compute the integral
    double prob = 0;

    // loop on the value of l
    for (int i = 0; i < l_size; i++) {

        // initalise the returned array to 0
        r_l[i] = 0;

        // loop to compute the integral
        for (int j = 0; j < n; j++){
            prob = exp(-(x[j] - w_r*w_r)*(x[j] - w_r*w_r) / 2.0 / w_r / w_r) / sqrt(2*M_PI) / w_r / w_r * w[j];
            r_l[i] = r_l[i] +  prob * s_l[i * n + j] * mu_I(x[j]);
        }
    }

    // ---------------------
    // freeing memory
    free_double_ptr(zeta);
    free_double_ptr(s_l);

    return r_l;
}




