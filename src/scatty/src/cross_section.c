#include "cross_section.h"
#include <stdio.h>


/* ----------------
compute the Coulomb phase shift for l and zeta = alpha / v_rel

delta_l = Ln[Gamma(1+l+I*zeta)/Gamma(1+l-I*zeta)]/(2*I) 
we implement the numerical scheme derived in arXiv:0912.3189 
---------------- */
double coulomb_phase_shift(int l, double zeta)
{

    double delta_0 = - zeta * GAMMA_E;
    
    double new_term = 0;
    int m = 0;

    do
    {
        new_term = ( atan(zeta/(m+1.0)) - zeta / (m+1.0) );
        delta_0 = delta_0 - new_term;
        m = m+1;

    } while (m < 1000 && fabs(new_term/delta_0) > 1e-4);
    
    if (l == 0)
        return delta_0;

    double delta_l = delta_0;

    for (m=1; m <=l; m++ )
        delta_l = delta_l + atan(zeta/(1.0*m));


    return delta_l;
}


double* coulomb_phase_shift_arr(int l, double* zeta_arr, int size)
{
    double* delta_ls = (double*)malloc(size * sizeof(double));
    if (!delta_ls) return NULL; // Handle memory allocation failure

    for (int i = 0; i < size; i++) {
        delta_ls[i] = coulomb_phase_shift(l, zeta_arr[i]);
    }

    return delta_ls;
}


void free_memory(double* ptr) {
    
    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
}
