import numpy as np

def periodic(w, S):
    """Applies periodic boundary conditions to
    postpoints

    inputs:
    w -- (ndarray, ndim=arbitrary) array whose indices are to be
         kept in range 0, 1, ... , S-1 per periodic BCs
    S -- (int) supremum value in the modular arithematic operation
         for example, S = z.N restricts indices between 0, 1, ... , z.N

    outputs:
    Array with periodic BCs enforced
    """

    return np.mod(w, S)

def periodic_old(z, i = None, Uf = None):
    """Applies periodic boundary conditions to
    postpoints

    inputs:
    z -- (instance) phase space variable
    i -- (int)
    Uf -- (ndarray, ndim=1) normalized fluxes
          for z.MCs

    outputs:
    k1, k2 -- (int) postpoint indices per
              periodic BCs
    """

    if i is None and Uf is None:
        # whole vector z.postpoints passed
        z.postpoints = np.mod(z.postpoints, z.N)
        return z.postpoints

    else:
        if Uf[i] > 0:
            k1 = z.postpoints[i]
            k2 = np.mod(z.postpoints[i] + 1, z.N)
        elif Uf[i] < 0:
            k1 = z.postpoints[i]
            k2 = np.mod(z.postpoints[i] - 1, z.N)
        else:
            # dummy entries, will store in summary
            # array izf as float('nan')
            k1 = None
            k2 = None
        return k1, k2
