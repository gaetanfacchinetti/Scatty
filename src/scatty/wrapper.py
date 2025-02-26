import os
import ctypes

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
LIB_PATH = os.path.join(BASE_DIR, 'src/libscatty.so')

LIB = ctypes.CDLL(LIB_PATH)

def wrap(func_name, args, res) -> None:

    func = getattr(LIB, func_name)
    
    func.argtypes = args
    func.restype = res

wrap('coulomb_phase_shift', [ctypes.c_int, ctypes.c_double], ctypes.c_double)
