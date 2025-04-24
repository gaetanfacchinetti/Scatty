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


#ifndef RATES_H
#define RATES_H

#include "integration.h"
#include "cross_section.h"

/**
 * @brief Weight of the energy transfer rate
 *
 * @param y  The angular momentum quantum number (int32_t*).
 * @return   ((y-1) + exp(-2*y)*(y+1)) / (y*y) (double).
 */
double mu_I(double y);


/**
 * @brief Partial wave contribution to the normalized energy transfer rate, r_\chi_l, 
 *        computed for the Coulomb phase shifts
 *
 * r_chi_l = 1/sqrt(2pi)/w_r^2 * int_0^\infty dy exp(-(y-w_r^2)^2/2/w_r^2) s_l(w_r zeta_r / y) mu(y)
 *
 * @param l     The angular momentum quantum number (int).
 * @param zeta  The zeta parameter (dimensionless) (double).
 * @return      Values of r_chi_l (double).
 */
double* r_chi_coulomb_arr(int32_t* l, double zeta_r, double w_r, int32_t l_size, double* x, double* w, int32_t n);


/**
 * @brief Partial wave contribution to the normalized energy transfer rate, r_\chi_l, 
 *        computed for the Coulomb phase shifts
 *
 * r_chi_l = 1/sqrt(2pi)/w_r^2 * int_0^\infty dy exp(-(y-w_r^2)^2/2/w_r^2) s_l(w_r zeta_r / y) mu(y)
 *
 * @param l     The angular momentum quantum number (int).
 * @param zeta  The zeta parameter (dimensionless) (double).
 * @return      Values of r_chi_l (double).
 */
double* r_chi_coulomb_grid(int32_t *l, double* zeta_r, double* w_r, int l_size, int32_t zeta_size, int32_t w_size, int32_t n);


double* r_chi_coulomb_ur_grid(int32_t*l, double* w_r, int32_t l_size, int32_t w_size);
double* r_chi_coulomb_mc_grid(int32_t *l, double* zeta_r, double* w_r, int32_t l_size, int32_t zeta_size, int32_t w_size, int32_t n);


#endif // RATES_H