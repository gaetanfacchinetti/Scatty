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

#include "tools.h"

void free_double_ptr(double* ptr) {
    
    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
}

void free_int_ptr(int* ptr) {
    
    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
}