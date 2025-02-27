import numpy as np
import ctypes


from .wrapper import LIB


def coulomb_phase_shift(l: int | np.ndarray, zeta: float | np.ndarray) -> np.ndarray:


    if isinstance(zeta, float):
        return LIB.coulomb_phase_shift(l, zeta)
    
    zeta_size = len(zeta)
    
    # Convert NumPy array to ctypes array
    zeta_ctypes = zeta.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    if isinstance(l, int):

        # Call the C function
        result_ptr = LIB.coulomb_phase_shift_arr(l, zeta_ctypes, zeta_size)

        # Convert C array back to NumPy array
        result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=(zeta_size,)))
   
    else:

        l_size = len(l)

        if np.any(np.diff(l) != 1):
            raise ValueError("Arrays of l must be given in ascending order and contiguous")
 
        # Convert NumPy array to ctypes array
        l_ctypes = np.array(l, dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    
        # Call the C function
        result_ptr = LIB.coulomb_phase_shift_grid(l_ctypes, zeta_ctypes, l_size, zeta_size)

        # Convert C array back to NumPy array
        result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=(l_size, zeta_size,)))


    # free the memory to avoid leaks
    LIB.free_double_ptr(result_ptr)

    return result_array