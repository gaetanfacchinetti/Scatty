import os
import ctypes

import numpy as np

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
LIB_PATH = os.path.join(BASE_DIR, 'src/libscatty.so')

LIB = ctypes.CDLL(LIB_PATH)

def wrap_func(func_name, args, res) -> None:

    func = getattr(LIB, func_name)
    
    func.argtypes = args
    func.restype = res


def wrap_args(*args):
    
    collapsed = []
    sizes = [None] * len(args)
    args_ctype = [None] * len(args)

    # loop over all arguments
    for i, arg in enumerate(args):

        if isinstance(arg, float | int):
            collapsed.append(i)
            arg = np.array([arg])

        if arg.dtype == np.int64:
            arg = np.array(arg, dtype=np.int32)

        if arg.dtype == np.float64:
            args_ctype[i] = arg.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        if arg.dtype == np.int32:
            args_ctype[i] = arg.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))

        sizes[i] = len(arg)

    return  tuple(sizes), tuple(collapsed), tuple(args_ctype)



wrap_func('free_double_ptr', [ctypes.POINTER(ctypes.c_double)], None)

wrap_func('coulomb_phase_shift', [ctypes.c_int32, ctypes.c_double], ctypes.c_double)
wrap_func('coulomb_phase_shift_grid', [ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_double), ctypes.c_int32, ctypes.c_int32], ctypes.POINTER(ctypes.c_double))

wrap_func('n_coulomb_transfer_cross_section_grid', [ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_double), ctypes.c_int32, ctypes.c_int], ctypes.POINTER(ctypes.c_double))
wrap_func('n_coulomb_ur_transfer_cross_section_arr', [ctypes.POINTER(ctypes.c_int32), ctypes.c_int32], ctypes.POINTER(ctypes.c_double))

#wrap('r_chi_coulomb_arr', [ctypes.POINTER(ctypes.c_int), ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int], ctypes.POINTER(ctypes.c_double))
wrap_func('r_chi_coulomb_grid', [ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int32,  ctypes.c_int32,  ctypes.c_int32, ctypes.c_int32], ctypes.POINTER(ctypes.c_double))
wrap_func('r_chi_coulomb_ur_grid', [ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_double), ctypes.c_int32, ctypes.c_int32], ctypes.POINTER(ctypes.c_double))


wrap_func('r_chi_coulomb_mc_grid', [ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int32,  ctypes.c_int32,  ctypes.c_int32, ctypes.c_int32], ctypes.POINTER(ctypes.c_double))

wrap_func('test_integral', [], ctypes.c_double)