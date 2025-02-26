import numpy as np
import ctypes


from .wrapper import LIB


def coulomb_phase_shift(l:int, zeta: float | np.ndarray):

    if isinstance(zeta, float):
        return LIB.coulomb_phase_shift(l, zeta)
    
    size = len(zeta)
    
    # Convert NumPy array to ctypes array
    zeta_ctypes = zeta.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    # Call the C function
    result_ptr = LIB.coulomb_phase_shift_arr(l, zeta_ctypes, size)

    # Convert C array back to NumPy array
    result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=(size,)))

    # free the memory to avoid leaks
    LIB.free_memory(result_ptr)

    return result_array