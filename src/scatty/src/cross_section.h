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

#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <math.h>
#include <string.h>

#include "tools.h"
#include "integration.h"
#include "cephes/complex.h"


// Euler-Mascheroni constant
#define GAMMA_E  0.57721566490153286

/**
 * @brief Phase shift of the Coulomb potential (for l and zeta = alpha / v_rel) 
 *
 * This function calculates the Coulomb phase shift using a given integer l
 * and a floating-point zeta parameter.
 * delta_l = Ln[Gamma(1+l+I*zeta)/Gamma(1+l-I*zeta)]/(2*I) 
 * Uses the cephes library https://www.netlib.org/cephes/
 * or the numerical scheme derived in arXiv:0912.3189 
 *
 * @param l     The angular momentum quantum number (int32_t).
 * @param zeta  The zeta parameter = alpha/v (dimensionless) (double).
 * @return      The Coulomb phase shift in radians (double).
 */
double coulomb_phase_shift(int32_t l, double zeta);


/**
 * @brief Phase shift of the Coulomb potential (for l and zeta = alpha / v_rel) 
 *
 * This function calculates the Coulomb phase shift using a given integer l
 * and a floating-point zeta parameter.
 * delta_l = Ln[Gamma(1+l+I*zeta)/Gamma(1+l-I*zeta)]/(2*I) 
 * Uses the cephes library https://www.netlib.org/cephes/
 * or the numerical scheme derived in arXiv:0912.3189 
 *
 * @param l       The angular momentum quantum number (int32_t*).
 * @param zeta    The zeta parameter = alpha/v (dimensionless) (double).
 * @param l_size  Size of the l array (int32_t).
 * @return        The Coulomb phase shift in radians (double).
 */
double* coulomb_phase_shift_arr(int32_t* l, double zeta, int32_t l_size);


/**
 * @brief Phase shift of the Coulomb potential (for l and zeta = alpha / v_rel) 
 *        where l and zeta are arrays (l being an ascending contiguous list of integers)
 *
 * This function calculates the Coulomb phase shift using a given integer l
 * and a floating-point zeta parameter.
 * delta_l = Ln[Gamma(1+l+I*zeta)/Gamma(1+l-I*zeta)]/(2*I) 
 * Uses the cephes library https://www.netlib.org/cephes/
 * or the numerical scheme derived in arXiv:0912.3189 
 *
 * @param l          The angular momentum quantum number (int32_t*).
 * @param zeta       The zeta parameter = alpha/v (dimensionless) (double*).
 * @param l_size     Size of the l array (int32_t).
 * @param zeta_size  Size of the zeta array (int32_t).
 * @return           The Coulomb phase shift in radians (double*).
 */
double* coulomb_phase_shift_grid(int32_t* l, double* zeta, int32_t l_size, int32_t zeta_size);


/**
 * @brief function of the phase shift appearing in the transfer cross-section
 *
 * t_l = sin(delta_l) * ((1+l) * sin(delta_l) + l * sin(delta_l - 2*delta_{l-1}))
 *
 * @param l          Angular momentum quantum number (int32_t).
 * @param delta_l    Phase shift for partial wave l (double).
 * @param delta_lm1  Phase shift for partial wave l-1 (double).
 * @return           Values of t_l (double).
 */
double transfer_factor(int32_t l, double delta_l, double delta_lm1);


/**
 * @brief Partial wave contribution to the normalized transfer cross-section, s_l, 
 *        computed for the Coulomb phase shifts
 *
 * By definition sigma_tansfer = 4\pi * zeta_4/(m_\chi^2 \alpha^2) sum_{l=0} s_l
 * (equivalent to sigma_tansfer = 4\pi alpha^2 /(m_\chi^2 v_rel^4) sum_{l=0} s_l).
 * s_l = sin(delta_l) * ((1+l) * sin(delta_l) + l * sin(delta_l - 2*delta_{l-1})) / zeta^2
 *
 * @param l     The angular momentum quantum number (int).
 * @param zeta  The zeta parameter (dimensionless) (double).
 * @return      Values of S_l (double).
 */
double n_coulomb_transfer_cross_section(int32_t l, double zeta);


/**
 * @brief Partial wave contribution to the normalized transfer cross-section, s_l, 
 *        computed for the Coulomb phase shifts
 *
 * By definition sigma_tansfer = 4\pi * zeta_4/(m_\chi^2 \alpha^2) sum_{l=0} s_l
 * (equivalent to sigma_tansfer = 4\pi alpha^2 /(m_\chi^2 v_rel^4) sum_{l=0} s_l).
 * s_l = sin(delta_l) * ((1+l) * sin(delta_l) + l * sin(delta_l - 2*delta_{l-1})) / zeta^2
 *
 * @param l         The angular momentum quantum number (int32_t*).
 * @param zeta      The zeta parameter (dimensionless) (double).
 * @param l_size    Size of the l array (int32_t).
 * @return          Values of S_l (double).
 */
double* n_coulomb_transfer_cross_section_arr(int32_t* l, double zeta, int32_t l_size);


/**
 * @brief Partial wave contribution to the normalized transfer cross-section, s_l, 
 *        computed for the Coulomb phase shifts
 *
 * By definition sigma_tansfer = 4\pi * zeta_4/(m_\chi^2 \alpha^2) sum_{l=0} s_l
 * (equivalent to sigma_tansfer = 4\pi alpha^2 /(m_\chi^2 v_rel^4) sum_{l=0} s_l).
 * s_l = sin(delta_l) * ((1+l) * sin(delta_l) + l * sin(delta_l - 2*delta_{l-1})) / zeta^2
 *
 * @param l          The angular momentum quantum number (int32_t*).
 * @param zeta       The zeta parameter (dimensionless) (double).
 * @param l_size     Size of the l array (int32_t).
 * @param zeta_size  Size of the zeta array (int32_t).
 * @return           Values of S_l (double).
 */
double* n_coulomb_transfer_cross_section_grid(int32_t *l, double* zeta, int32_t l_size, int32_t zeta_size);


/**
 * @brief Partial wave contribution to the normalized transfer cross-section, s_l, 
 *        computed for the Coulomb phase shifts with small coupling alpha << v
 *
 * By definition sigma_tansfer = 4\pi * zeta_4/(m_\chi^2 \alpha^2) sum_{l=0} s_l
 * (equivalent to sigma_tansfer = 4\pi alpha^2 /(m_\chi^2 v_rel^4) sum_{l=0} s_l).
 * s_l = sin(delta_l) * ((1+l) * sin(delta_l) + l * sin(delta_l - 2*delta_{l-1})) / zeta^2
 *
 * @param l         The angular momentum quantum number (int32_t*).
 * @param l_size    Size of the l array (int32_t).
 * @return          Values of S_l (double).
 */
double* n_coulomb_ur_transfer_cross_section_arr(int32_t* l, int32_t l_size);



#endif // CROSS_SECTION_H