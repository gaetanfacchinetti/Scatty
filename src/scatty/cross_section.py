import numpy as np
import ctypes


from .wrapper import LIB


def coulomb_phase_shift(l: int | np.ndarray, zeta: float | np.ndarray) -> np.ndarray:

    not_arr = np.array([], dtype=int)

    if isinstance(l, int):
        not_arr = np.append(not_arr, 0)
        l = np.array([l], dtype=np.int32)

    if isinstance(zeta, float):
        not_arr = np.append(not_arr, 1)
        zeta = np.array([zeta])

   
    zeta_size = len(zeta)
    l_size    = len(l)

    # convert NumPy array to ctypes array
    zeta_ctypes = zeta.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    l_ctypes    = np.array(l, dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    # call the C function and convert C array back to NumPy array
    result_ptr   = result_ptr = LIB.coulomb_phase_shift_grid(l_ctypes, zeta_ctypes, l_size, zeta_size)
    result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=(l_size, zeta_size,)))
    
    # free the result pointer
    LIB.free_double_ptr(result_ptr)

    result_array = result_array.squeeze(axis=tuple(not_arr))

    return result_array



def n_coulomb_transfer_cross_section(l: int | np.ndarray, zeta: float | np.ndarray) -> np.ndarray:
    
    not_arr = np.array([], dtype=int)

    if isinstance(l, int):
        not_arr = np.append(not_arr, 0)
        l = np.array([l], dtype=np.int32)

    if isinstance(zeta, float):
        not_arr = np.append(not_arr, 1)
        zeta = np.array([zeta])

   
    zeta_size = len(zeta)
    l_size    = len(l)

    # convert NumPy array to ctypes array
    zeta_ctypes = zeta.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    l_ctypes    = np.array(l, dtype=np.int32).ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    # call the C function and convert C array back to NumPy array
    result_ptr   = LIB.n_coulomb_transfer_cross_section_grid(l_ctypes, zeta_ctypes, l_size, zeta_size)
    result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=(l_size, zeta_size,)))
    
    # free the result pointer
    LIB.free_double_ptr(result_ptr)

    #result_array = result_array.squeeze(axis=tuple(not_arr))

    return result_array


def test_integral():
    
    res = LIB.test_integral()
    return res
    