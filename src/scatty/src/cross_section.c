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

#include "cross_section.h"

/** -----------------------------------------------
 * PHASE SHIFTS
 */

double coulomb_phase_shift(int32_t l, double zeta)
{

    // convention
    if (l == -1) return 0.0;

    // use cephes library for precision and rapidity
    double res = creal(  clog( cgamma(1.0 + l + zeta*I) / cgamma(1.0 + l - zeta*I) ) / (2.0*I)  );
    
    // if cephes fails use a more time consuming method
    if (!(isnan(res))){
        return res;
    } else {

        // first compute delta for l = 0
        double delta_0 = - zeta * GAMMA_E;
        
        double new_term = 0;
        int32_t m = 0;

        do {
            new_term = ( atan(zeta/(m+1.0)) - zeta / (m+1.0) );
            delta_0 = delta_0 - new_term;
            m = m+1;

        } while (m < 100000 && fabs(new_term/delta_0) > 1e-12);
        
        if (l == 0)
            return delta_0;

        // if l > 0 need to add more terms
        double delta_l = delta_0;

        for (m=1; m <=l; m++){
            delta_l = delta_l + zeta/(1.0*m);
        }

        return delta_l;
    }

    return NAN;

}


double* coulomb_phase_shift_arr(int32_t* l, double zeta, int32_t l_size)
{
    double* delta_l = (double*)malloc(l_size * sizeof(double));
    if (!delta_l) return NULL; // handle memory allocation failure

    int i_init = 0;

    // special case where the series starts with lm1
    if (l[0] == -1) {
        i_init = 1;
        delta_l[0] =  zeta > 0 ? 0.0 : NAN;
    }
        
    // first get the values for the minimum l
    // assumes the values of l are given in ascending order
    delta_l[i_init] = zeta > 0 ? coulomb_phase_shift(l[i_init], zeta) : NAN;

  
    // then keeps increasing for highest values of l
    // values of l must be given as a contiguous list
    for (int32_t i = i_init + 1; i < l_size; i++) {
        delta_l[i] = zeta > 0 ? delta_l[i-1] +  atan(zeta/(1.0*l[i])) : NAN;
    }

    return delta_l;
}


double* coulomb_phase_shift_grid(int32_t* l, double* zeta, int32_t l_size, int32_t zeta_size)
{
    double* delta_l = (double*)malloc(l_size * zeta_size * sizeof(double));
    if (!delta_l) return NULL; // handle memory allocation failure

    int i_init = 0;

    // special case where the series starts with lm1
    if (l[0] == -1) {
        i_init = 1;
        for (int32_t j = 0; j < zeta_size; j++) {
            delta_l[j] =  zeta[j] > 0 ? 0.0 : NAN;
        }
    }
        
    // first get the values for the minimum l
    // assumes the values of l are given in ascending order
    for (int32_t j = 0; j < zeta_size; j++) {
        delta_l[i_init * zeta_size + j] = zeta[j] > 0 ? coulomb_phase_shift(l[i_init], zeta[j]) : NAN;
    }
  
    // then keeps increasing for highest values of l
    // values of l must be given as a contiguous list
    for (int32_t i = i_init + 1; i < l_size; i++) {
        for (int32_t j = 0; j < zeta_size; j++){
            delta_l[i * zeta_size + j] = zeta[j] > 0 ? delta_l[(i-1)*zeta_size + j] +  atan(zeta[j]/(1.0*l[i])) : NAN;
        }
    }

    return delta_l;
}

/** -----------------------------------------------
 * TRANSFER CROSS SECTIONS
 */


double transfer_factor(int32_t l, double delta_l, double delta_lm1) 
{
    return sin(delta_l) * ( (1+l) * sin(delta_l) + l * sin(delta_l - 2*delta_lm1) );
}


double n_coulomb_transfer_cross_section(int32_t l, double zeta)
{
    // only compute the phase shift once and get the adjacent values fast
    double delta_lm1 = l > 0 ? coulomb_phase_shift(l-1, zeta) : 0.0;
    double delta_l   = l > 0 ? delta_lm1 + atan(zeta/(1.0*l)) : coulomb_phase_shift(0, zeta);
    
    return transfer_factor(l, delta_l, delta_lm1) / zeta / zeta;
}

double* n_coulomb_transfer_cross_section_arr(int32_t* l, double zeta, int32_t l_size)
{
    // first get the values of the phase shifts
    int32_t* lm1 = (int*)malloc((l_size+1) * sizeof(int32_t));
    if (!lm1) return NULL;

    lm1[0] = l[0]-1;
    memmove(lm1 + 1, l, l_size * sizeof(int32_t));

    double* delta_lm1 = coulomb_phase_shift_arr(lm1, zeta, l_size+1);
    
    // then compute the array of normalised l-wave cross-section 
    double* s_l = (double*)malloc(l_size * sizeof(double));
    if (!s_l) return NULL;

    for (int32_t i = 0; i < l_size; i++) {
        s_l[i] = transfer_factor(l[i], delta_lm1[(i+1)], delta_lm1[i]) / zeta / zeta;
    }

    // free the arrays of ls and deltas
    free_int_ptr(lm1);
    free_double_ptr(delta_lm1);

    return s_l;
}


double* n_coulomb_transfer_cross_section_grid(int32_t* l, double* zeta, int32_t l_size, int32_t zeta_size)
{
    // first get the values of the phase shifts
    // as we need l-1 they are all "shifted" to the right
    // so need an array that has 1 additional memory block
    int32_t* lm1 = (int*)malloc((l_size+1) * sizeof(int32_t));
    if (!lm1) return NULL;

    // shifting the values
    lm1[0] = l[0]-1;
    memmove(lm1 + 1, l, l_size * sizeof(int32_t));

    double* delta_lm1 = coulomb_phase_shift_grid(lm1, zeta, l_size+1, zeta_size);
    
    // then compute the array of normalised l-wave cross-section 
    double* s_l = (double*)malloc(l_size * zeta_size * sizeof(double));
    if (!s_l) return NULL;

    for (int32_t i = 0; i < l_size; i++) {
        for (int32_t j = 0; j < zeta_size; j++){
            s_l[i * zeta_size + j] = transfer_factor(l[i], delta_lm1[(i+1) * zeta_size + j], delta_lm1[i * zeta_size + j]) / zeta[j] / zeta[j];
        }
    }

    // free the arrays of ls and deltas
    free_int_ptr(lm1);
    free_double_ptr(delta_lm1);

    return s_l;
}


// s_l for zeta = 0
// l array should be in increasing order (but not necessarily continous)
double* n_coulomb_ur_transfer_cross_section_arr(int32_t* l, int32_t l_size)
{
    double* s_l = (double*)malloc(l_size * sizeof(double));
    if (!s_l) return NULL;

    int m = l[l_size - 1], j = 0;

    if (l[0] == 0){
        s_l[j] = pow(GAMMA_E, 2);
        j = j+1;
    }

    // initialise the digamma function at 1 (l=0, l+1=1)
    double psi = - GAMMA_E;

    for (int32_t i = 1; i <= m; i++){
        // psi(1+l) = psi(l) + 1/i
        psi = psi + 1.0/i;

        if (i == l[j]){
            s_l[j] = psi*(psi+2.0);
            j = j+1;
        }
    } 

    return s_l;
}







