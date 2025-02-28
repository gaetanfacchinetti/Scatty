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