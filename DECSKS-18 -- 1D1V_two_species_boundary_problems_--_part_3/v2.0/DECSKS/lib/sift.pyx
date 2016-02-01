import numpy as np
cimport numpy as np

ctypedef np.float64_t DTYPE_t

def sift(np.ndarray[DTYPE_t, ndim=2] in, np.ndarray[DTYPE_t, ndim=2] vx):

    cdef int i
    cdef int j
    cdef int Nx = vx.shape[0]
    cdef int Nvx  = vx.shape[1]
    cdef list out = []

    for i in range(Nx):
        for j in range(Nvx):
            if vx[i,j] >= 0:
                out.append(in[i,j])

    return out
                
