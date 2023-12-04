#/usr/bin/env python3

import numpy as np
from thermodynamics_parameters import R, beta, omega0, alpha


def omega(d):
    """
    Loop entropy function under the exact powerlaw model.
    """
    
    return omega0*d**(-alpha)


def W(i,V10):
    """
    Scaler product of exact loop entropy.
    """
    
    j = np.arange(2, i)
    omegas_array = omega(2*(i + 1 - j))
    
    
    return np.array(list(V10.values()))[1:-1]@omegas_array