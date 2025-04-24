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


#include "rates.h"

/** -----------------------------------------------
 * ENERGY TRANSFER RATES
 */


double mu_I(double y)
{
    return ((y-1) + exp(-2*y)*(y+1)) / (y*y);
}

double* r_chi_coulomb_arr(int32_t* l, double zeta_r, double w_r, int32_t l_size, double* x, double* w, int32_t n)
{
    // --------------------------
    // set the weights and roots
    double* x_loc = (double*)malloc(n * sizeof(double));
    double* w_loc = (double*)malloc(n * sizeof(double));

    if (!x_loc) return NULL;
    if (!w_loc) return NULL;

    memcpy(x_loc, x, n * sizeof(double));
    memcpy(w_loc, w, n * sizeof(double));

    int err = set_root_and_weights_scale(n, fmax(w_r*(-5.0 + w_r), 0.0), w_r*(5.0 + w_r), x_loc, w_loc, false);
    if (err < 0) return NULL;

    // --------------------------
    // prepare the integration

    // get the grid of values on which to evaluate the cross-section
    double* zeta = (double*)malloc(n*sizeof(double));
    if (!zeta) return NULL;

    for (int32_t j=0; j< n; j++) 
        zeta[j] = w_r * zeta_r / x_loc[j];

    // get the values of s_l we need to integrate over
    double* s_l = n_coulomb_transfer_cross_section_grid(l, zeta, l_size, n);

    // preparing the returned array
    double* r_l = (double*)calloc(l_size,  sizeof(double));
    if (!r_l) return NULL;

    // ----------------------------
    // compute the integral
    double prob = 0.0, mu = 0.0;

    // loop to compute the integral
    for (int32_t j = 0; j < n; j++){

        prob = exp(- pow(x_loc[j]/w_r - w_r, 2) / 2.0) / sqrt(2*M_PI) / w_r / w_r  * w_loc[j];
        mu = mu_I(x_loc[j]);
        
        // loop on the value of l
        for (int32_t i = 0; i < l_size; i++) {
            r_l[i] += prob * s_l[i * n + j] * mu;
        }

    }

    // ---------------------
    // free memory
    free_double_ptr(zeta);
    free_double_ptr(s_l);
    free_double_ptr(x_loc);
    free_double_ptr(w_loc);

    return r_l;
}


double* r_chi_coulomb_grid(int32_t *l, double* zeta_r, double* w_r, int32_t l_size, int32_t zeta_size, int32_t w_size, int32_t n)
{

    double* x = (double*)malloc(n * sizeof(double));
    double* w = (double*)malloc(n * sizeof(double));

    if (!x) return NULL;
    if (!w) return NULL;

    int err = set_gauss_legendre_points_and_weights(n, x, w);
    if (err < 0) return NULL;

    double* r_chi = (double*)malloc(l_size * zeta_size * w_size * sizeof(double));
    if (!r_chi) return NULL; // handle memory allocation failure

    double* r_chi_arr;
    
    // NEED TO REARRANGE THE ORDER FOR FASTER EVALUATION
    for (int32_t j = 0; j < zeta_size; j++){

        for (int32_t k = 0; k < w_size; k++){

            // call the function array
            r_chi_arr = r_chi_coulomb_arr(l, zeta_r[j], w_r[k], l_size, x, w, n);
            
            for (int32_t i = 0; i < l_size; i++) {
                r_chi[(i * zeta_size + j) * w_size + k] = r_chi_arr[i];

            }
        }
    }

    // ---------------------
    // free memory
    free_double_ptr(r_chi_arr);
    free_double_ptr(x);
    free_double_ptr(w);

    return r_chi;
    
}


double* r_chi_coulomb_ur_grid(int32_t*l, double* w_r, int32_t l_size, int32_t w_size)
{
    // get the values of s_l we need to integrate over
    double* s_l = n_coulomb_ur_transfer_cross_section_arr(l, l_size);

    double* r_l = (double*)malloc(l_size * w_size * sizeof(double));
    if (!r_l) return NULL; // handle memory allocation failure

    double part1 = 0, part2 = 0;

    for (int32_t j = 0; j < w_size; j++){
        
        part1 = erf(w_r[j]/sqrt(2.0));
        part2 = - sqrt(2.0/M_PI) * w_r[j] * exp(-w_r[j]*w_r[j]/2.0);
        
        for (int32_t i = 0; i < l_size; i++) {
            r_l[i * w_size + j] = s_l[i] * (part1 + part2)/ pow(w_r[j], 3);
            //r_l[i * w_size + j] = s_l[i] ;//* (part1 + part2)/ pow(w_r[j], 3);
        }
    }

    // ---------------------
    // free memory
    free_double_ptr(s_l);

    return r_l;
}



double* r_chi_coulomb_mc_grid(int32_t *l, double* zeta_r, double* w_r, int32_t l_size, int32_t zeta_size, int w_size, int n)
{

    double* x = (double*)malloc(n * sizeof(double));
    if (!x) return NULL;

    draw_gauss_monte_carlo_points(n, x, 0.0, 1.0);

    double* r_chi = (double*)malloc(l_size * zeta_size * w_size * sizeof(double));
    if (!r_chi) return NULL; // handle memory allocation failure

    double prob = 0.0, y = 0.0;
    double* s_l = NULL;
    
    for (int32_t j = 0; j < zeta_size; j++){
        for (int32_t k = 0; k < w_size; k++){
            for (int32_t m = 0; m < n; m++){ 

                y = w_r[k] * (w_r[k] + x[m]);
                s_l = n_coulomb_transfer_cross_section_arr(l, w_r[k] * zeta_r[j] / y, l_size);
                    
                for (int32_t i = 0; i < l_size; i++){
                    // loop on the value of l
                    r_chi[(j * w_size + k) * l_size + i] +=  (y > 0 ? s_l[i] * mu_I(y) / w_r[k]  / n : 0.0);

                }
            }
        }
    }

    // ---------------------
    // free memory
    free_double_ptr(x);
    free_double_ptr(s_l);

    return r_chi;
    
}