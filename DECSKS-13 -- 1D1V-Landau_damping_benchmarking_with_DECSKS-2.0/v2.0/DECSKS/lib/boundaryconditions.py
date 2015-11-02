import numpy as np

def periodic(w, Nw):
    """Applies periodic boundary conditions to
    postpoints

    inputs:
    w -- (ndarray, ndim=arbitrary) array whose indices are to be
         kept in range 0, 1, ... , S-1 per periodic BCs
    S -- (int) supremum value in the modular arithematic operation
         for example, S = z.N restricts indices between 0, 1, ... , z.N

    outputs:
    Array with periodic BCs enforced according to bound on attribute4 S.N
    """

    return np.mod(w, Nw)
