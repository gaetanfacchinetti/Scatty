#include "cross_section.h"

double dummy(double x)
{
    return x;
}

/* compute the Coulomb phase shift for l and zeta = alpha / v_rel
delta_l = Ln[Gamma(1+l+I*zeta)/Gamma(1+l-I*zeta)]/(2*I) */
double c_coulomb_phase_shift(int l, double zeta)
{
    // Calculate the gamma function
    gsl_sf_result norm, arg1, arg2; 

    gsl_sf_lngamma_complex_e(1.0 + l, zeta, &norm, &arg1);
    gsl_sf_lngamma_complex_e(1.0 + l, -zeta, &norm, &arg2);
    
    /* 
    we know that the ratio of the two Gamma functions is of norm 1
    therefore, only the difference of the arguments matters
    each of them can be in ]-pi, pi], therefore the difference is in ]-2*pi, 2*pi]
    we bring it back to ]-pi, pi] with a translation
    */
    double res = arg1.val - arg2.val;

    res = (res > M_PI) ? res - 2*M_PI : res;
    res = (res <= -M_PI) ? res + 2*M_PI : res;

    // divide by two by definition and get a result in ]-pi/2, pi/2]
    res = res / 2.0;
        
    return res;
}
