#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <math.h>
#include <string.h>

#include "tools.h"
#include "integration.h"

#define GAMMA_E  0.57721566490153286;

/**
 * @brief Phase shift of the Coulomb potential (for l and zeta = alpha / v_rel) 
 *
 * This function calculates the Coulomb phase shift using a given integer l
 * and a floating-point zeta parameter.
 * delta_l = Ln[Gamma(1+l+I*zeta)/Gamma(1+l-I*zeta)]/(2*I) 
 * We implement the numerical scheme derived in arXiv:0912.3189 
 *
 * @param l     The angular momentum quantum number (int).
 * @param zeta  The zeta parameter (dimensionless) (double).
 * @return      The Coulomb phase shift in radians (double).
 */
double coulomb_phase_shift(int l, double zeta);

/**
 * @brief Phase shift of the Coulomb potential (for l and zeta = alpha / v_rel) 
 *        where zeta is an array
 *
 * This function calculates the Coulomb phase shift using a given integer l
 * and a floating-point zeta parameter.
 * delta_l = Ln[Gamma(1+l+I*zeta)/Gamma(1+l-I*zeta)]/(2*I) 
 * we implement the numerical scheme derived in arXiv:0912.3189 
 *
 * @param l     The angular momentum quantum number (int).
 * @param zeta  The zeta parameter (dimensionless) (double*).
 * @return      The Coulomb phase shift in radians (double*).
 */
double* coulomb_phase_shift_arr(int l, double* zeta, int zeta_size);

/**
 * @brief Phase shift of the Coulomb potential (for l and zeta = alpha / v_rel) 
 *        where l and zeta are arrays (l being an ascending contiguous list of integers)
 *
 * This function calculates the Coulomb phase shift using a given list of integers l
 * and a given list of zeta.
 * delta_l = Ln[Gamma(1+l+I*zeta)/Gamma(1+l-I*zeta)]/(2*I) 
 * we implement the numerical scheme derived in arXiv:0912.3189 
 *
 * @param l     The angular momentum quantum number (int *).
 * @param zeta  The zeta parameter (dimensionless) (double*).
 * @return      The Coulomb phase shift in radians (double*).
 */
double* coulomb_phase_shift_grid(int* l, double* zeta, int l_size, int zeta_size);


/**
 * @brief Partial wave contribution to the normalized transfer cross-section, S_l, 
 *        computed for the Coulomb phase shifts
 *
 * By definition sigma_tansfer = 4\pi * zeta_4/(m_\chi^2 \alpha^2) sum_{l=0} S_l
 * (equivalent to sigma_tansfer = 4\pi alpha^2 /(m_\chi^2 v_rel^4) sum_{l=0} S_l).
 * S_l = sin(delta_l) * ((1+l) * sin(delta_l) + l * sin(delta_l - 2*delta_{l-1})) / zeta^2
 *
 * @param l     The angular momentum quantum number (int).
 * @param zeta  The zeta parameter (dimensionless) (double).
 * @return      Values of S_l (double).
 */
double n_coulomb_transfer_cross_section(int l, double zeta);


double* n_coulomb_transfer_cross_section_grid(int *l, double* zeta, int l_size, int zeta_size);

double transfer_factor(int l, double delta_l, double delta_lm1);

double mu_I(double y);
double r_chi_coulomb(int l, double zeta_r, double w_r);

#endif // CROSS_SECTION_H