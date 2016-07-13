import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
DTYPEINT = np.int64

ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t DTYPEINT_t

@cython.boundscheck(False)
def create_postpoint_objects(np.ndarray[DTYPE_t, ndim=1] S,
                             np.ndarray[DTYPEINT_t, ndim=2] GCFL_int,
                             np.ndarray[DTYPE_t, ndim=2] GCFL_frac,
                             int num_stages,
                             int Nx,
                             int Nvx,
                             DTYPE_t xwidth,
                             DTYPE_t vxwidth,
                             sim_params
                             ):
    """
    create objects f_S and N_S that are the objects that capture
    the effect of a bulk source flowing in from the right

    inputs:
    S -- (ndarray, ndim=1) source function data, contains only data for vx < 0
    GCFL_int -- (ndarray, ndim=2) indexed [s,j], int(CFL) for stage s at vx j
    GCFL_frac -- (ndarray, ndim=2) indexed [s,j], frac(CFL) for stage s at vx j
    num _stages -- (int) sim_params['splitting']['number_of_stages']['a'],
                   total number of stages

    Nx  -- (int) total number of x gridpoints
    Nvx -- (int) total number of vx gridpoints
    sim_params -- (dict) simulation paramters

    outputs:
    f_S -- (ndarray, ndim=3), shape = (stages+1, x.N, vx.N) density array that is
           the contribution of the source particles to every [i,j] for
           a velocity dependent distribution function

    N_S -- (ndarray, ndim = 1), shape = (stages+1,) the number of source
            particles that reach the lower boundary

    see notebook s28 for details
    """
    cdef int S_Nvx = len(S)
    cdef int s, i, j
    cdef int k1, k2

    # initialize an array is indexed as (s,i,j) where s labels the stage (>= 1)
    cdef np.ndarray[DTYPE_t, ndim=3] f_S = np.zeros((num_stages + 1, Nx, Nvx))

    # initialize a 1D array indexed by s, the stage (s >= 1; s = 0 entry is empty)
    cdef np.ndarray[DTYPE_t, ndim=1] N_S = np.zeros(num_stages + 1)

    for s in range(1, num_stages + 1):
        for j in range(S_Nvx):

            NG = np.abs(GCFL_int[s,j]) + 1 # each j requires NG ghost cells
            g_crit = np.abs(GCFL_int[s,j]) - Nx # the critical ghost index g (in configuration, i = Nx + g) for a given j where g < g_crit reaches the wall, g >= g_crit maps on-grid

            for g in range(NG): # indexes position of ghost cell: i = Nx + g
                # note these cells are not stored, but their effect has been
                # worked out, see notebook s28

                if 0 <= g < g_crit: # both proportions k1 and k2 reach boundary
                    N_S[s] +=  S[j] * (1 + GCFL_frac[s,j])
                    N_S[s] += -S[j] *  GCFL_frac[s,j]

                elif g == g_crit: # edge case, ghost prepoints reach i in [-1,0]
                                  # lower boundary at i = -1/2

                    if -GCFL_frac[s,j] > 1/2.: # both k1 and k2 reach boundary
                        # physically k1 = -1/2, k2 = -1
                        N_S[s] +=  S[j] * (1 + GCFL_frac[s,j])
                        N_S[s] += -S[j] *  GCFL_frac[s,j]

                    elif -GCFL_frac[s,j] < 1/2.: # k2 reaches boundary, k1 on grid
                        # physically k1 = 0, k2 = -1/2

                        f_S[s,k1,j] += S[j] * (1 + 2 * GCFL_frac[s,j])
                        N_S[s] += -2 * S[j] *  GCFL_frac[s,j]

                elif np.max((g_crit + 1, 0)) <= g < NG-1:

                    # both proportions map on-grid to points (g + GCFL.int[s,j]) and (g + GCFL.int[s,j] - 1)
                    k1 = g + GCFL_int[s,j]
                    k2 = k1 - 1

                    f_S[s,k1,j] +=  S[j] * (1 + GCFL_frac[s,j])
                    f_S[s,k2,j] += -S[j] *  GCFL_frac[s,j]

                elif g == NG-1:
                    # only k2 proportion maps on-grid
                    # physically k2 = Nx-1, k1 = Nx (off-grid, ghost region)

                    k2 = Nx-1
                    f_S[s,k2,j] += -S[j] *  GCFL_frac[s,j]

                else:
                    # this case (g > NG-1) never happens, loop is for 0 <= g <= NG-1, included for aesthetics
                    pass

    # up to this point total number that has reached the wall has been stored in
    # a 1D vector N_S[s] for each stage, where N_S[s] = np.sum(chi[k]*S[s,j])
    # need to integrate over space:
    #
    #         $$N_S = \int dx \int dv_x S(x_g -> x[k] <= x_wall)
    #
    #               \simeq \sum_k \sum_g \sum_j \Delta x \Delta v_x chi_k S
    #
    # factor in multiplicative prefactors to all entry N_S[s]:
    N_S *= xwidth
    N_S *= vxwidth

    return f_S, N_S
