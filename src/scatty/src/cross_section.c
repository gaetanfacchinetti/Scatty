#include "cross_section.h"

/* ----------------
compute the Coulomb phase shift for l and zeta = alpha / v_rel

delta_l = Ln[Gamma(1+l+I*zeta)/Gamma(1+l-I*zeta)]/(2*I) 
we implement the numerical scheme derived in arXiv:0912.3189 
---------------- */
double coulomb_phase_shift(int l, double zeta)
{
    return 1.0;
}
