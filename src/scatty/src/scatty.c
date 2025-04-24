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

int main()
{

    // testing the implementation of the phase shifts
    double dl_0 = coulomb_phase_shift(0, 0.1);
    double dl_1 = coulomb_phase_shift(0, 10.0);
    double dl_2 = coulomb_phase_shift(1, 0.1);

    printf("Test of phase shifts: %f, %f, %f \n", dl_0, dl_1, dl_2);


    return 0;
}