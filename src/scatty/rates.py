import numpy as np

from .wrapper import LIB, wrap_args


# normalised Coulomb transfer cross-section
def r_chi_coulomb(l: int | np.ndarray, 
                  zeta_r: float | np.ndarray, 
                  w_r: float | np.ndarray,
                  n: int = 5000, 
                  method: str = 'MC') -> np.ndarray:
    """
    Energy transfer rate induced by scattering in coulomb potential

    Parameters
    ----------
    l : int | np.ndarray
        index of the partial wave expansion
    zeta_r : float | np.ndarray
        dimensionless ration of the coupling strength 
        over the velocity of the scattered particle alpha / v / u_r
        with u_r the relative velocity dispersion of the particles
    w_r : float | np.ndarray
        ratio of the bulk relative velocity over the relative velocity dispersion
    n : int, optional
        number of integration point. Default is 5000.
    method : str, optional
        Integration method, can be MC for Monte-Carlo or GL for Gauss-Legendre
        Default is 'MC'
        
    Returns
    -------
    np.ndarray

    Raises
    ------
    ValueError
        when method is set to somethig else than MC or GL
    """

    # prepare the arguments to wrap
    sizes, collapsed, args = wrap_args(l, zeta_r, w_r)

    if method == 'GL':

        # call the C function and convert C array back to NumPy array
        result_ptr   = LIB.r_chi_coulomb_grid(*args, *sizes, n)
        result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=sizes))

        # free the C pointer
        LIB.free_double_ptr(result_ptr)
        return result_array.squeeze(collapsed)

    if method == 'MC':

        # prepare the arguments to wrap
        sizes, collapsed, args = wrap_args(l, zeta_r, w_r)

        # call the C function and convert C array back to NumPy array
        # need to be careful here as the table is indexes from zeta, w, l
        # and NOT l, zeta, w (need to perform a permutation)
        result_ptr   = LIB.r_chi_coulomb_mc_grid(*args, *sizes, n)
        result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=(sizes[1], sizes[2], sizes[0],)))

        # free the C pointer 
        LIB.free_double_ptr(result_ptr)

        return result_array.transpose((2, 0, 1)).squeeze(axis=tuple(collapsed))
    
    else:
        raise ValueError("Need to specify method of integration 'GL' or 'MC'.")



# normalised Coulomb transfer cross-section
def r_chi_coulomb_ur(l: int | np.ndarray, 
                     w_r: float | np.ndarray) -> np.ndarray:
    """
    Energy transfer rate induced by scattering in coulomb potential
    for small coupling constant zeta, zeta -> 0

    Parameters
    ----------
    l : int | np.ndarray
        index of the partial wave expansion
    w_r : float | np.ndarray
        ratio of the bulk relative velocity over the relative velocity dispersion

    Returns
    -------
    np.ndarray
    """
    
    # prepare the arguments to wrap
    sizes, collapsed, args = wrap_args(l, w_r)

    # call the C function and convert C array back to NumPy array
    result_ptr   = LIB.r_chi_coulomb_ur_grid(*args, *sizes)
    result_array = np.copy(np.ctypeslib.as_array(result_ptr, shape=sizes))

    # free the C pointer
    LIB.free_double_ptr(result_ptr)

    return result_array.squeeze(collapsed)


