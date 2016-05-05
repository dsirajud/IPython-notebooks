import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
DTYPEINT = np.int64

ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t DTYPEINT_t


def nearest_gridpoint_assignment(np.ndarray[DTYPE_t, ndim=2] f_old,  np.ndarray[DTYPE_t, ndim=2] Uf_old, np.ndarray[DTYPEINT_t, ndim=2] map1, np.ndarray[DTYPEINT_t, ndim=2] map2, int dim1, int dim2):

    cdef np.ndarray[DTYPE_t, ndim=2] f_new = np.zeros([dim1, dim2], dtype = DTYPE)
    cdef int i
    cdef int j

    for i in range(dim1):
        for j in range(dim2):
            if Uf_old[i,j] < 0:
                f_new[map1[i,j], map2[i,j]] += (f_old[i,j] + Uf_old[i,j])
            else: # Uf_old[i,j] is nonnegative
                f_new[map1[i,j], map2[i,j]] += (f_old[i,j] - Uf_old[i,j])
    return f_new

def contiguous_gridpoint_assignment(np.ndarray[DTYPE_t, ndim=2] f_old,  np.ndarray[DTYPE_t, ndim=2] Uf_old, np.ndarray[DTYPEINT_t, ndim=2] map1, np.ndarray[DTYPEINT_t, ndim=2] map2, int dim1, int dim2):

    cdef np.ndarray[DTYPE_t, ndim=2] f_new = np.zeros([dim1, dim2], dtype = DTYPE)
    cdef int i
    cdef int j

    for i in range(dim1):
        for j in range(dim2):
            if Uf_old[i,j] < 0:
                f_new[map1[i,j], map2[i,j]] -= Uf_old[i,j]
            else: # Uf_old[i,j] is nonnegative
                f_new[map1[i,j], map2[i,j]] += Uf_old[i,j]
    return f_new

def assignment(np.ndarray[DTYPE_t, ndim=2] f_old,  np.ndarray[DTYPEINT_t, ndim=2] map1, np.ndarray[DTYPEINT_t, ndim=2] map2, int dim1, int dim2):

    cdef np.ndarray[DTYPE_t, ndim=2] f_new = np.zeros([dim1, dim2], dtype = DTYPE)
    cdef int i
    cdef int j

    for i in range(dim1):
        for j in range(dim2):
            f_new[map1[i,j], map2[i,j]] += f_old[i,j]
    return f_new

def sift(np.ndarray[DTYPE_t, ndim=2] Uf):
    """
    inputs:
    Uf -- (ndarray, ndim=2, dtype = float64) flux like object

    output:
    Uf_nonneg, Uf_neg -- (ndarrays, ndim=2, dtype = float64) containers holding
                    flux entires from Uf according to

                    pos -- for each f[i,j] such that its Uf[i,j] >= 0
                    neg -- for each f[i,j] such that its Uf[i,j] < 0

    """
    cdef int i, j
    cdef int dim1 = Uf.shape[0]
    cdef int dim2 = Uf.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] Uf_nonneg = np.zeros((dim1, dim2))
    cdef np.ndarray[DTYPE_t, ndim=2] Uf_neg = np.zeros((dim1, dim2))

    for j in range(dim2):
        for i in range(dim1):
            if Uf[i,j] >= 0:
                for i in range(dim1):
                    Uf_nonneg[i,j] = Uf[i,j]
            else:
                for i in range(dim1):
                    Uf_neg[i,j] = Uf[i,j]

    return Uf_nonneg, Uf_neg

def sift_by_column(np.ndarray[DTYPE_t, ndim=2] Uf):
    """
    Assumes Uf[i,j] has the same sign for all i at a given column j.
    this routine was coded for boundary conditions that did not change
    the flux array Uf entries (i.e. any boundary condition except
    a symmetric boundary condition). For a symmetry boundary, should
    use the general sift method above.

    inputs:
    Uf -- (ndarray, ndim=2, dtype = float64) flux like object

    output:
    Uf_nonneg, Uf_neg -- (ndarrays, ndim=2, dtype = float64) containers holding
                    flux entires from Uf according to

                    pos -- for each f[i,j] such that its Uf[i,j] >= 0
                    neg -- for each f[i,j] such that its Uf[i,j] < 0

    """
    cdef int i, j
    cdef int dim1 = Uf.shape[0]
    cdef int dim2 = Uf.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] Uf_nonneg = np.zeros((dim1, dim2))
    cdef np.ndarray[DTYPE_t, ndim=2] Uf_neg = np.zeros((dim1, dim2))

    for j in range(dim2):
        # CFL[0,j] = CFL[1,j] = ... = const, we choose to examine (0,j) arbitrarily
        if Uf[0,j] >= 0:
            for i in range(dim1):
                Uf_nonneg[i,j] = Uf[i,j]
        else:
            for i in range(dim1):
                Uf_neg[i,j] = Uf[i,j]

    return Uf_nonneg, Uf_neg
