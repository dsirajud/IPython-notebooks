import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
DTYPEINT = np.int64

ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t DTYPEINT_t


def assignment(np.ndarray[DTYPE_t, ndim=2] f_old,  np.ndarray[DTYPEINT_t, ndim=2] map1, np.ndarray[DTYPEINT_t, ndim=2] map2, int dim1, int dim2):

    cdef np.ndarray[DTYPE_t, ndim=2] f_new = np.zeros([dim1, dim2], dtype = DTYPE)
    cdef int i
    cdef int j

    for i in range(dim1):
        for j in range(dim2):
            f_new[map1[i,j], map2[i,j]] += f_old[i,j]
    return f_new

def sift(np.ndarray[DTYPE_t, ndim=2] f, np.ndarray[DTYPE_t, ndim=2] CFL):
    """
    inputs:
    f -- (ndarray, ndim=2, dtype = float64) density like object
    vz -- (ndarray, ndim=2, dtype = float64) generalized velocity grid

    output:
    f_nonneg, f_neg -- (ndarrays, ndim=2, dtype = float64) containers holding
                    density like entires from f according to

                    pos -- for each f[i,j] such that its vz[i,j] >= 0
                    neg -- for each f[i,j] such that its vz[i,j] < 0

    """
    cdef int i, j
    cdef int dim1 = f.shape[0]
    cdef int dim2 = f.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] f_nonneg = np.zeros((dim1, dim2))
    cdef np.ndarray[DTYPE_t, ndim=2] f_neg = np.zeros((dim1, dim2))

    for i in range(dim1):
        for j in range(dim2):
            if CFL[i,j] >= 0:
                f_nonneg[i,j] = f[i,j]
            else:
                f_neg[i,j] = f[i,j]

    return f_nonneg, f_neg
