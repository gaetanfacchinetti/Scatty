import os
import ctypes

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
LIB_PATH = os.path.join(BASE_DIR, 'src/libscatty.so')

LIB = ctypes.CDLL(LIB_PATH)

def wrap(func_name, args, res) -> None:

    func = getattr(LIB, func_name)
    
    func.argtypes = args
    func.restype = res


wrap('free_double_ptr', [ctypes.POINTER(ctypes.c_double)], None)

wrap('coulomb_phase_shift', [ctypes.c_int, ctypes.c_double], ctypes.c_double)
wrap('coulomb_phase_shift_arr', [ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_int], ctypes.POINTER(ctypes.c_double))
wrap('coulomb_phase_shift_grid', [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int], ctypes.POINTER(ctypes.c_double))

wrap('n_coulomb_transfer_cross_section_grid', [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int], ctypes.POINTER(ctypes.c_double))
wrap('test_integral', [], ctypes.c_double)