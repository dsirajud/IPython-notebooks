import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
DTYPEINT = np.int64

ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t DTYPEINT_t

def flux(np.ndarray[DTYPE_t, ndim=2] c,
         np.ndarray[DTYPE_t, ndim=3] d,
         np.ndarray[DTYPEINT_t, ndim=2] vz_true_prepointmesh)

    cdef int dim0 = d.shape[0]
    cdef int dim1 = d.shape[1]
    cdef int dim2 = d.shape[2]
    cdef np.ndarray[DTYPE_t, ndim=2] Uf = np.zeros((dim1, dim2))
    cdef int i, j, q
    for j in range(dim2):
        for i in range(dim1):
            if j == vz_true_prepointmesh[i,j]:
                for q in range(dim0):
                    Uf[i,j] += c[q,j]*d[q,i,j]
                else:
                    Uf[i,j] += c[q,vz_true_prepointmesh[i,j]] \
                      * (-1) ** q * d[q,i,j]

    return Uf
