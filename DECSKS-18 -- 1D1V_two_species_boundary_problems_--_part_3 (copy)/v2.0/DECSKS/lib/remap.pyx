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

    for j in range(dim2):
        # CFL[0,j] = CFL[1,j] = ... = const, we choose to examine (0,j) arbitrarily
        if CFL[0,j] >= 0:
            for i in range(dim1):
                f_nonneg[i,j] = f[i,j]
        else:
            for i in range(dim1):
                f_neg[i,j] = f[i,j]

    return f_nonneg, f_neg

def symmetric(np.ndarray[DTYPE_t, ndim=2] f,
              np.ndarray[DTYPE_t, ndim=2] Uf,
              np.ndarray[DTYPEINT_t, ndim=2] map1,
              np.ndarray[DTYPEINT_t, ndim=2] map2,
              int Nx, int Nvx):

    cdef np.ndarray[DTYPE_t, ndim=2] f_new = np.zeros([Nx, Nvx], dtype = DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] Uf_star = np.zeros([Nx, Nvx], dtype = DTYPE)
    cdef int i
    cdef int j

    # TODO map1 is not k1, and map2 is not k2, rather map1 is z.postpointmesh[k,:,:] (for a passed value of k = 0, 1) and vz.prepointmesh is map2
    # TODO require passing the 3D matrices in order to access all mappings here.

    for i in range(Nx):
        for j in range(Nvx):

            if map1[i,j] <= (Nx - 1):
                # apply right boundary condition for partners

            elif map1[i,j] == -(Nx - 2):
                k1_star = -map1[i,j]
                k2_star = -map2[i,j] - 1
                j_star = Nvx - 1 - j
                # compute Uf_star

                f_new[k1_star, j_star] += f[i,j] - Uf_star[i,j]

                # apply right BC for partner at k2_star

            elif -(Nx - 3) <= map1[i,j] <= -1:
                k1_star = -map1[i,j]
                k2_star = -map2[i,j] - 1
                j_star = Nvx - 1 - j

                # compute Uf_star
                f_new[k1_star, j_star] += f[i,j] - Uf_star[i,j]
                f_new[k2_star, j_star] += Uf_star[i,j]


            elif map1[i,j] == 0:
                k1_star = -map1[i,j]
                k2_star = -map2[i,j] - 1
                j_star = Nvx - 1 - j
                # compute Uf_star

                f_new[ map1[i,j], map2[i,j] ] += f[i,j] + Uf[i,j]
                f_new[ k1_star, j_star] += f[i,j] - Uf_star[i,j]
                f_new[ k2_star, j_star] += Uf_star[i,j]

            elif map1[i,j] == 1:
                k1_star = -map1[i,j]
                k2_star = -map2[i,j] - 1
                j_star = Nvx - 1 - j
                # compute Uf_star

                f_new[ map1[i,j], map2[i,j] ] += f[i,j] + Uf[i,j]
                f_new[ k1_star, j_star] += f[i,j] - Uf_star[i,j]
                f_new[ k2_star, j_star] += Uf_star[i,j]


    return f_new
