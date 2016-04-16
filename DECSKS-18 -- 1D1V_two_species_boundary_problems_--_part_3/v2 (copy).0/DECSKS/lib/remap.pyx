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

def sift_flux(np.ndarray[DTYPE_t, ndim=2] Uf):
    """
    Examines a flux-like array Uf, and returns two arrays, Uf_nonneg
    and Uf_neg. These two arrays which are identically sized with respect
    to Uf and are so-named due to Uf_nonneg containing the entries of
    Uf which are nonnegative, and Uf_neg cdontaining those which
    are negative. Note that this method expects a flux-like object (Uf),
    not a density-like object (f). Recall that,

        Uf = CFL.frac * f + (Higher order corrections)

    since f is strictly nonnegative everywhere in phase space, the
    flux Uf already factors in CFL.frac as the zeroeth order term,
    and the numerical limiter on the high order corrections ensure that Uf
    shares the same sign as the zeroeth order term, hence Uf has the same
    sign as CFL.frac.

    This routine differs from remap.sift_old in three ways:

        (1) it sifts a flux Uf, rather than a density f
        (2) it checks the sign on the input paramter Uf, hence
            checking the sign on CFL.frac is not needed as discussed
            above.
        (3) it checks the sign of every (i,j), not just the first row
            entry in a column.

    check on every (i,j), whereas remap.sift_old only checks the top
    entry for each j (i.e. is used for situations where each column has
    one CFL number per column j that corresponds to all i therein). Hence,
    this routine (remap.sift) is general, whereas the other (remap.sift_old)
    is only valid for some cases, but not all (so far, the only boundary
    condition that precludes using remap.sift_old is the symmetry/reflecting
    boundary condition)
    """
    cdef int i, j
    cdef int dim1 = Uf.shape[0]
    cdef int dim2 = Uf.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=2] Uf_nonneg = np.zeros((dim1, dim2))
    cdef np.ndarray[DTYPE_t, ndim=2] Uf_neg = np.zeros((dim1, dim2))

    for j in range(dim2):
        for i in range(dim1):
            if Uf[i,j] >= 0:
                Uf_nonneg[i,j] = Uf[i,j]
            else:
                Uf_neg[i,j] = Uf[i,j]

    return Uf_nonneg, Uf_neg

def sift_old(np.ndarray[DTYPE_t, ndim=2] f, np.ndarray[DTYPE_t, ndim=2] CFL):
    """
    Examines a density-like array f, and returns two arrays, f_nonneg
    and f_neg. These two arrays which are identically sized with respect
    to f and are so-named due to f_nonneg containing the entries of
    f which correspond to CFL.numbers[i,j] with nonnegative sign, and f_neg
    containing those which have negative sign.

    inputs:

    f -- (ndarray, ndim = 2) density-like or flux-like array
    CFL -- (ndarray, ndim = 2) 2D array of CFL numbers, identical dimenions as f

    outputs:

    f_nonneg -- (ndarray, ndim = 2) density-like array with only those entries
                from f[i,j] which are associated with sign(CFL.numbers[i,j]) >= 0

    f_neg -- (ndarray, ndim = 2) density-like array with only those entries
                from f[i,j] which are associated with sign(CFL.numbers[i,j]) < 0


    NOTE:

    This routine is of limited applicability. It presumes that the sign of
    the normalized velocity (CFL number) only varies in the column direction,
    i.e. that each row within a column is described by the same CFL number
    which is the physical velocity normalized by the grid velocity:

       CFL number = v / v_grid = v / (x.width / t.width) = v * t.width / x.width

    This will be true for all cases EXCEPT a symmetric/reflection boundary condition
    which reverses the direction of any such density packet which reaches the
    so-described boundary in a time step. Since not all rows i in the
    prepoint set {(i,j)} at a given j will reach the boundary (where
    (i,j) refer to configuration and velocity prepoints, respectively), only
    some rows i, in CFL[i,j] will have their direction reversed. Hence, such a
    boundary condition breaks the constancy among rows for a given column.

    This routine exploits the constancy of CFL sign among the rows for each column
    to gain efficiency. We only check the sign of the first CFL value in the column
    and presume all rows share the same sign as the first row. Again, this will
    be true for all BCs except a symmetry/reflecting condition. If used, it reduces
    the number of conditional checks on the sign of every (i,j) (requiring Nx*Nvx
    checks) to conditional checks on the sign of one j for all such i (requires only
    Nvx checks).

    The above routine (remap.sift) checks every (i,j) and hence is general.
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
