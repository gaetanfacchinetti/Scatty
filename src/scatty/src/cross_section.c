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
    double* s_l     = (double*)malloc(l_size * zeta_size * sizeof(double));
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


double r_chi_coulomb_integrand(double y, int l, double zeta_r, double w_r)
{
    return 0.0;
}

double mu_I(double y)
{
    return ((y-1) + exp(-2*y)*(y+1)) / (y*y);
}