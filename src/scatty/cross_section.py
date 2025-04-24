import numpy as np

from .wrapper import LIB, wrap_args


def coulomb_phase_shift(l: int | np.ndarray, 
                        zeta: float | np.ndarray) -> np.ndarray:
    """
    Phase shift induced by scattering in coulomb potential

    Parameters
    ----------
    l : int | np.ndarray
        index of the partial wave expansion
    zeta : float | np.ndarray
        dimensionless ration of the coupling strength 
        over the velocity of the scattered particle alpha / v

    Returns
    -------
    np.ndarray
    """

    # prepare the arguments to wrap
    sizes, collapsed, args = wrap_args(l, zeta)

    # call the C function and convert C array back to NumPy array
    result_ptr   = LIB.coulomb_phase_shift_grid(*args, *sizes)
    result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=sizes))
    
    # free the result pointer
    LIB.free_double_ptr(result_ptr)

    return result_array.squeeze(axis=collapsed)



# normalised Coulomb transfer cross-section
def n_coulomb_transfer_cross_section(l: int | np.ndarray, 
                                     zeta: float | np.ndarray) -> np.ndarray:
    """
    Transfer cross-sectrion induced by scattering in coulomb potential

    Parameters
    ----------
    l : int | np.ndarray
        index of the partial wave expansion
    zeta : float | np.ndarray
        dimensionless ratio of the coupling strength 
        over the velocity of the scattered particle alpha / v

    Returns
    -------
    np.ndarray
    """

    # prepare the arguments to wrap
    sizes, collapsed, args = wrap_args(l, zeta)

    # call the C function and convert C array back to NumPy array
    result_ptr   = LIB.n_coulomb_transfer_cross_section_grid(*args, *sizes)
    result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=sizes))
    
    # free the result pointer
    LIB.free_double_ptr(result_ptr)

    return result_array.squeeze(axis=collapsed)




# normalised Coulomb transfer cross-section (for zeta = 0, alpha /  v -> 0)
def n_coulomb_ur_transfer_cross_section(l: int | np.ndarray) -> np.ndarray:
    """
    Transfer cross-sectrion induced by scattering in coulomb potential
    with small coupling constant zeta 

    Parameters
    ----------
    l : int | np.ndarray
        index of the partial wave expansion
    zeta : float | np.ndarray
        dimensionless ratio of the coupling strength 
        over the velocity of the scattered particle alpha / v

    Returns
    -------
    np.ndarray
    """
    
    # prepare the arguments to wrap
    sizes, collapsed, args = wrap_args(l)
    
    # call the C function and convert C array back to NumPy array
    result_ptr   = LIB.n_coulomb_ur_transfer_cross_section_arr(*args, *sizes)
    result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=sizes))
    
    # free the result pointer
    LIB.free_double_ptr(result_ptr)

    return result_array.squeeze(axis=collapsed)




def test_integral():
    
    res = LIB.test_integral()
    return res
    