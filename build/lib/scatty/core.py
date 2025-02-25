import ctypes

def test():
    
    # Load the shared library
    scatty_lib = ctypes.CDLL('../../src_c/libscatty.so')

    # Define the function prototypes if necessary
    scatty_lib.dummy.argtypes = [ctypes.c_double]
    scatty_lib.dummy.restype = ctypes.c_double

    # Call the function
    result = scatty_lib.dummy(1.0)
    print(result)