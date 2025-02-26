#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <math.h>
#include <stdlib.h>

#define GAMMA_E  0.57721566490153286;

void free_memory(double* ptr);
double coulomb_phase_shift(int l, double zeta);
double* coulomb_phase_shift_arr(int l, double* zeta_arr, int size);

#endif // CROSS_SECTION_H