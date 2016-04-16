import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
DTYPEINT = np.int64

ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t DTYPEINT_t

# PYTHON METHODS
def periodic(f_old,
             Uf,
             z,
             vz,
             sim_params,
             charge,
             k = 0
             ):
    """Applies periodic boundary conditions to
    postpointmesh

    inputs:
    f_old -- (ndarray, ndim=2) density array
    z -- (instance) phase space variable being evolved

        z.postpointmesh -- (ndarray, ndim=3),
                           shape = (2, x.N, vx.N)

    outputs:
    f_old -- (ndarray, ndim=2) Array with both periodic
             BCs being enforce
    z    -- (instance) phase sapce variable being evolved with
             updated attribute z.postpointmesh

    f_old, Uf returned for symmetry with nonperiodic routine below
    """
    z.postpointmesh[k,:,:] = np.mod(z.postpointmesh[k,:,:], z.N)
    vz.postpointmesh[k,:,:] = vz.prepointmesh

    return f_old, Uf

def nonperiodic(f_old,
                Uf,
                z,
                vz,
                sim_params,
                charge,
                k = 0
                ):
    """orchestrates applying nonperiodic boundary conditions
    to the array w with total active grid points Nw. Nonperiodic
    boundary conditions require handling both left and right
    boundaries

    inputs:
    f_old -- (ndarray, ndim=2) density array
    z -- (instance) phase space variable being evolved

    outputs:
    f_old -- (ndarray, ndim=2) density with both left and right
             nonperiodic BCs enforced
    Uf -- (ndarray, ndim=2) high order fluxes with both left and right
             nonperiodic BCs enforced


    z returned (no changes) for symmetry with periodic routine above
    """
    # lower boundary
    f_old, Uf = eval(sim_params['BC'][z.str]['lower'] +
                           '_lower_boundary')(f_old,
                                              Uf,
                                              z.postpointmesh[k,:,:],
                                              vz.prepointmesh,
                                              z.N, vz.N, k, charge,
                                              sim_params,
                                              z, vz)

    # upper boundary
    f_old, Uf = eval(sim_params['BC'][z.str]['upper'] +
                           '_upper_boundary')(f_old,
                                              Uf,
                                              z.postpointmesh[k,:,:],
                                              vz.prepointmesh,
                                              z.N, vz.N, k, charge,
                                              sim_params,
                                              z, vz)


    vz.postpointmesh[k,:,:] = vz.prepointmesh

    return f_old, Uf

# CYTHON METHODS
@cython.boundscheck(False)
def absorbing_lower_boundary(np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              int Nz, int Nvz, int k, int charge,
                              sim_params,
                              z, vz):

    # vars here are typed as C data types to minimize python interaction
    cdef int i, j
    for i in range(Nz):
        for j in range(Nvz):
            if zpostpointmesh[i,j] <= 0:
                f_old[i,j] = 0
                Uf_old[i,j] = 0
                zpostpointmesh[i,j] = 0 # set postpoint at the absorber

    # permanently copy to instance attribute
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old

@cython.boundscheck(False)
def absorbing_upper_boundary(np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              int Nz, int Nvz, int k, int charge,
                              sim_params,
                              z, vz):

    # vars here are typed as C data types to minimize python interaction
    cdef int i, j
    for i in range(Nz):
        for j in range(Nvz):
            if zpostpointmesh[i,j] >= Nz - 1:
                f_old[i,j] = 0
                Uf_old[i,j] = 0

                zpostpointmesh[i,j] = Nz - 1 # set postpoint at the absorber

    # permanently copy to instance attribute
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old

@cython.boundscheck(False)
def collector_lower_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointvaluemesh,
        int Nz, int Nvz, int k, int charge,
        sim_params,
        z, vz
        ):
    # vars here are typed as C data types to minimize python interaction

    # for any such density packet whose postpoint predicts exceeding the lower boundary of
    # the domain, zero out (absorb) the density, and add contribution to total charge
    # density sigma at that boundary
    cdef int i, j, j_inflection_point
    cdef DTYPE_t sigma_nk = 0
    cdef DTYPE_t vzwidth = vz.width

    j_inflection = sim_params['vx_inflection_gridpoint']
    if k == 0:
        for i in range(Nz):
            for j in range(0,j_inflection):
                if zpostpointmesh[i,j] <= 0:
                    sigma_nk += -vzprepointvaluemesh[i,j] * (f_old[i,j] + Uf_old[i,j])
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = Nz - 1  # set postpoint at the absorber

    elif k == 1:
        for i in range(Nz):
            for j in range(0,j_inflection):
                if zpostpointmesh[i,j] <= 0:
                    sigma_nk += vzprepointvaluemesh[i,j] * Uf_old[i,j]
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = Nz - 1  # set postpoint at the absorber


    sigma_nk *= charge
    sigma_nk *= vzwidth

    # update cumulative charge density
    sim_params['sigma'][z.str]['lower'] += sigma_nk
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old

@cython.boundscheck(False)
def collector_upper_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointvaluemesh,
        int Nz, int Nvz, int k, int charge,
        sim_params,
        z, vz
        ):
    # vars here are typed as C data types to minimize python interaction

    # for any such density packet whose postpoint predicts exceeding the lower boundary of
    # the domain, zero out (absorb) the density, and add contribution to total charge
    # density sigma at that boundary
    cdef int i, j, j_inflection_point
    cdef DTYPE_t sigma_nk = 0
    cdef DTYPE_t vzwidth = vz.width

    j_inflection = sim_params['vx_inflection_gridpoint']

    if k == 0:
        for i in range(Nz):
            for j in range(j_inflection,Nvz):
                if zpostpointmesh[i,j] >=  Nz - 1:
                    sigma_nk += vzprepointvaluemesh[i,j] * (f_old[i,j] - Uf_old[i,j])
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = 0  # set postpoint at the absorber

    elif k == 1:
        for i in range(Nz):
            for j in range(j_inflection,Nvz):
                if zpostpointmesh[i,j] >=  Nz - 1:
                    sigma_nk += vzprepointvaluemesh[i,j] * Uf_old[i,j]
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = 0  # set postpoint at the absorber


    sigma_nk *= charge
    sigma_nk *= vzwidth

    # update cumulative charge density
    sim_params['sigma'][z.str]['upper'] += sigma_nk
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old
