#include "cross_section.h"


int main()
{

    // testing the implementation of the phase shifts
    double dl_0 = coulomb_phase_shift(0, 0.1);
    double dl_1 = coulomb_phase_shift(0, 10.0);
    double dl_2 = coulomb_phase_shift(1, 0.1);

    printf("Test of phase shifts: %f, %f, %f \n", dl_0, dl_1, dl_2);


    return 0;
}