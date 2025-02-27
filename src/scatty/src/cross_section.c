#include "cross_section.h"

double coulomb_phase_shift(int l, double zeta)
{
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
    double* delta_ls = (double*)malloc(zeta_size * sizeof(double));
    if (!delta_ls) return NULL; // handle memory allocation failure

    // first get the values for the minimum l
    // assumes the values of l are given in ascending order
    for (int j = 0; j < zeta_size ; j++) {
        delta_ls[j] = coulomb_phase_shift(l, zeta[j]);
    }
            
    return delta_ls;
}


double* coulomb_phase_shift_grid(int* l, double* zeta, int l_size, int zeta_size)
{
    double* delta_ls = (double*)malloc(l_size * zeta_size * sizeof(double));
    if (!delta_ls) return NULL; // handle memory allocation failure

    // first get the values for the minimum l
    // assumes the values of l are given in ascending order
    for (int j = 0; j < zeta_size; j++) {
        delta_ls[j] = coulomb_phase_shift(l[0], zeta[j]);
    }
  
    // then keeps increasing for highest values of l
    // values of l must be given as a contiguous list
    for (int i = 1; i < l_size; i++) {
        for (int j = 0; j < zeta_size; j++){
            delta_ls[i * zeta_size + j] = delta_ls[(i-1)*zeta_size + j] +  atan(zeta[j]/(1.0*l[i]));
        }
    }

    return delta_ls;
}


double n_coulomb_transfer_cross_section(int l, double zeta)
{
    // only compute the phase shift once and get the adjacent values fast
    double delta_lm1 = l > 0 ? coulomb_phase_shift(l-1, zeta) : 0.0;
    double delta_l   = l > 0 ? delta_lm1 + atan(zeta/(1.0*l)) : coulomb_phase_shift(0, zeta);
    
    return sin(delta_l) * ( (1+l) * sin(delta_l) + l * sin(delta_l - 2*delta_lm1) ) / zeta / zeta;
}


double* n_coulomb_transfer_cross_section_grid(int *l, double* zeta, int l_size, int zeta_size)
{
    double* sl = (double*)malloc(l_size * zeta_size * sizeof(double));
    if (!sl) return NULL; // handle memory allocation failure

    double* delta_lm1 = (double*)malloc(l_size * zeta_size * sizeof(double));
    if (!delta_lm1) return NULL;

    if (l[0] > 0) {

        int *lm1 = (int *)malloc(l_size * sizeof(int));
        if (!lm1) return NULL;

        for (int i = 0; i < l_size; i++)
            lm1[i] = l[i] -1;
        
        delta_lm1 = coulomb_phase_shift_grid(lm1, zeta, l_size, zeta_size);
    }

    //

    return sl;

}


double mu_I(double y)
{
    return ((y-1) + exp(-2*y)*(y+1)) / (y*y);
}