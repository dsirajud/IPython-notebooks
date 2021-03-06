import DECSKS
import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
DTYPEINT = np.int64

ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t DTYPEINT_t

# PYTHON METHODS
def periodic(f_k,
             f_old,
             Uf,
             t,
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
    vz.postpointmesh[k,:,:] = vz.prepointmesh.copy()   # assign to a copy so that changes to postpointmesh do not affect prepointmesh

    return f_k, f_old, Uf, z

def nonperiodic(f_k,
                f_old,
                Uf,
                dt,
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

    # initialize vz.postpointmesh to be prepointmesh values. Some boundary routines may require modifying these, most will not
    vz.postpointmesh[k,:,:] = vz.prepointmesh.copy() # assign to a copy so that changes to postpointmesh do not affect prepointmesh

    # lower boundary
    f_k, f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['lower'])(
                                              f_k,     # only lower symmetry boundary would need access to this
                                              f_old,
                                              Uf,
                                              z.postpointmesh[k,:,:],
                                              vz.postpointmesh[k,:,:],
                                              vz.prepointmesh,
                                              vz.prepointvaluemesh,
                                              dt,
                                              z.N, vz.N, k, charge,
                                              sim_params,
                                              z, vz)

    # upper boundary
    f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['upper'])(
                                              f_old,
                                              Uf,
                                              z.postpointmesh[k,:,:],
                                              vz.postpointmesh[k,:,:],
                                              vz.prepointmesh,
                                              vz.prepointvaluemesh,
                                              dt,
                                              z.N, vz.N, k, charge,
                                              sim_params,
                                              z, vz)


    #    vz.postpointmesh[k,:,:] = vz.prepointmesh
    return f_k, f_old, Uf, z

# CYTHON METHODS
@cython.boundscheck(False)
def absorbing_lower_boundary(np.ndarray[DTYPE_t, ndim=2] f_k,
                              np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
                              float dt,
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

    return f_k, f_old, Uf_old

@cython.boundscheck(False)
def absorbing_upper_boundary(np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
                              float dt,
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
def cutoff_lower_boundary(np.ndarray[DTYPE_t, ndim=2] f_k,
                          np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
                              float dt,
                              int Nz, int Nvz, int k, int charge,
                              sim_params,
                              z, vz):

    # vars here are typed as C data types to minimize python interaction
    cdef int i, j
    for i in range(Nz):
        for j in range(Nvz):
            if zpostpointmesh[i,j] < 0:
                f_old[i,j] = 0
                Uf_old[i,j] = 0
                zpostpointmesh[i,j] = 0 # set postpoint whose density/flux is zero
                                        # to a dummy gridpoint, e.g. lower boundary

    # permanently copy to instance attribute
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_k, f_old, Uf_old

@cython.boundscheck(False)
def cutoff_upper_boundary(np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
                              float dt,
                              int Nz, int Nvz, int k, int charge,
                              sim_params,
                              z, vz):

    # vars here are typed as C data types to minimize python interaction
    cdef int i, j
    for i in range(Nz):
        for j in range(Nvz):
            if zpostpointmesh[i,j] > Nz - 1:
                f_old[i,j] = 0
                Uf_old[i,j] = 0

                zpostpointmesh[i,j] = Nz - 1 # set postpoint whose density/flux is zero
                                             # to a dummy gridpoint, e.g. upper boundary

    # permanently copy to instance attribute
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old

def collector_lower_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_k,
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
        float dt,
        int Nz, int Nvz, int k, int charge,
        sim_params,
        z, vz
        ):
    # vars here are typed as C data types to minimize python interaction

    # for any such density packet whose postpoint predicts exceeding the lower boundary of
    # the domain, zero out (absorb) the density, and add contribution to total charge
    # density sigma at that boundary
    cdef int i, j
    cdef DTYPE_t sigma_nk = 0
    cdef DTYPE_t vzwidth = vz.width


    if k == 0:
        for i in range(Nz):
            for j in range(Nvz):
                # while the flux limiter ensures the signs of vz.prepointvaluemesh[i,j] is consistent with
                # Uf_old[i,j], it is more natural here to ascertain the sign on vz.prepointvaluemesh[i,j]
                # as we would like to tack on contributions from fluxes
                #
                #    vz.prepoitnvaluemesh[i,j] * (f_old[i,j] - Uf_old[i,j])
                #
                #           or
                #
                #    vz.prepointvaluemesh[i,j]* Uf_old[i,j]
                #
                # note here we store objects like f_old and Uf_old, but the equivalent meaning of
                # (f - Uf) = (1 - U)f, i.e. some fraction of the distribution function f.
                # It is just an implementation decision to not compute U directly, but instead retain
                # the product Uf.
                if zpostpointmesh[i,j] <= 0 and vzprepointvaluemesh[i,j] < 0:
                    sigma_nk += (f_old[i,j] + Uf_old[i,j])
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = 0  # set postpoint at the absorber

    if k == 1:
        for i in range(Nz):
            for j in range(Nvz):
                if zpostpointmesh[i,j] <= 0 and vzprepointvaluemesh[i,j] < 0:
                    sigma_nk +=  (-Uf_old[i,j])
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = 0  # set postpoint at the absorber


    sigma_nk *= charge
    sigma_nk *= vzwidth
    sim_params['sigma'][z.str]['lower'] += sigma_nk

    z.postpointmesh[k,:,:] = zpostpointmesh


    return f_k, f_old, Uf_old

def collector_upper_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
        float dt,
        int Nz, int Nvz, int k, int charge,
        sim_params,
        z, vz
        ):
    # vars here are typed as C data types to minimize python interaction

    # for any such density packet whose postpoint predicts exceeding the lower boundary of
    # the domain, zero out (absorb) the density, and add contribution to total charge
    # density sigma at that boundary
    cdef int i, j
    cdef DTYPE_t sigma_nk = 0
    cdef DTYPE_t vzwidth = vz.width

    if k == 0:
        for i in range(Nz):
            for j in range(Nvz): # rightward velocities (incoming fluxes) at the right boundary
                if zpostpointmesh[i,j] >=  Nz - 1:
                    if vzprepointvaluemesh[i,j] > 0:
                        sigma_nk += (f_old[i,j] + Uf_old[i,j])
                        f_old[i,j] = 0
                        Uf_old[i,j] = 0
                        zpostpointmesh[i,j] = Nz - 2  # set postpoint at the absorber
                    else:
                        zpostpointmesh[i,j] = Nz - 2  # set postpoint at the absorber


    if k == 1:
        for i in range(Nz):
            for j in range(Nvz): # rightward velocities (incoming fluxes) at the right boundary
                if zpostpointmesh[i,j] >=  Nz - 1:
                    if vzprepointvaluemesh[i,j] > 0:
                        sigma_nk += -Uf_old[i,j]
                        Uf_old[i,j] = 0
                        zpostpointmesh[i,j] = Nz - 2  # set postpoint at the absorber
                    else:
                        zpostpointmesh[i,j] = Nz - 2  # set postpoint at the absorber

    # update cumulative charge density sigma = sum_k sum_j sum_n sigma_nk * dvx * dt
    # k = k1, k2, j = 0, ..., Nvx -1, n = 0, 1, 2, ... Nt

    sigma_nk *= charge
    sigma_nk *= vzwidth
    sim_params['sigma'][z.str]['upper'] += sigma_nk
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old

@cython.boundscheck(False)
def symmetric_lower_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_k,
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
        float dt,
        int Nz, int Nvz, int k,
        int charge,
        sim_params,
        z, vz
        ):

    """
    This routine assumes a time substep is positive as the Vlasov-Poisson system is irreversible
    when boundaries are present (i.e. the upper boundary cannot be periodic/open in this setup)

    this boundary condition must not inherently operate on a symmetric velocity grid,
    but it is more natural and consistent this way. Some rationale is given below:

    For example, if an exiting particle at the symmetry plane has
    a velocity V1 < 0, it is the physical case that there exists an entering partner particle from the
    complement domain with oppositely directly velocity -V1 > 0. If a symmetric velocity
    grid, i.e. [avx, bvx] = [-V, V], is used then it is always possible to factor in both
    exiting (V1 < 0) density with entering partners (-V1 > 0) provided that V1 is inside [-V,V].

    If we have an asymmetric case where [avx, bvx] where avx != -bvx, then in the handling
    of the symmetry plane, there would be cases where such an exiting particle having V1 < 0 would have
    its entering partner with velocity -V1 not be representable (i.e. if avx < V1 is in the velocity grid,
    but bvx < -V1). We permit options so that the user may elect to take the following alternative to
    a symmetric velocity grid, i.e.

    We could still argue that it would be appropriate to handle such situations where the entering partner
    velocity is off-grid (at too largely positive of a value for a given grid [avx, bvx])
    then to take the partner density to be zero and interpret the  situation as the "window" or
    control volume in phase space. Thus, any entering particle with off-grid velocities will not re-enter
    the domain as their phase space coordinates are not in the simulation domain.

    a CUTOFF upper boundary condition should be chosen for the distribution function on the associated
    velocity variable for any case. CAUTION: periodic BCs on velocity for the distribution function are not
    appropriate, as the reflection in velocity exiting grid point j -> vz.N - 1 - j will not successfully mirror
    to the exact grid point needed if vz.N != vz.Ngridpoints (as would be the case for periodic boundary conditions)
    """
    vzpostpointmesh = vzprepointmesh.copy() # by this point, vzpostpointmesh has not been assigned (is an array of zeroes), init to prepointmesh values
                                            # and replace exiting particles with their partner velocities if encountered when scanning of all [i,j]

    cdef int i, j
    if k == 0:
        for j in range(Nvz/2):
            for i in range(Nz):
                if zpostpointmesh[i,j] <= 0:
                    # prepare for partner remap in lib.convect_configuration.remap_step
                    zpostpointmesh[i,j] = -zpostpointmesh[i,j]
                    vzpostpointmesh[i,j] = Nvz - 1 - vzprepointmesh[i,j]
                    Uf_old[i,j] = -Uf_old[i,j]

    if k == 1:
        for j in range(Nvz / 2):
            for i in range(Nz):
                if zpostpointmesh[i,j] <= 0:
                    # prepare for partner remap in lib.convect_configuration.remap_step
                    zpostpointmesh[i,j] = -zpostpointmesh[i,j]
                    vzpostpointmesh[i,j] = Nvz - 1 - vzpostpointmesh[i,j]
                    Uf_old[i,j] = -Uf_old[i,j]

    z.postpointmesh[k,:,:] = zpostpointmesh
    vz.postpointmesh[k,:,:] = vzpostpointmesh

    return f_k, f_old, Uf_old
