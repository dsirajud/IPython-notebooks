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

    NOTE: f_k is the distribution function contribution to all postpoints k.
    Only a symmetry boundary condition requires access to this object as it needs
    to make sure particles whose postpoints are exactly on the symmetry plane
    are factored in before replacing them in the reused f_old and Uf_old objects
    with their entering partners which are subsequently remapped as usual. This
    special consideration is described in github.com/dsirajud/../DECSKS-18 part 3 notebook.
    """

    # this orchestrator is called if periodic boundary conditions are specified on the distribution
    # function if a velocity variable or if periodic boundary conditions are specfied on BOTH
    # the distribution function for the configuration variable and its corresponding boundary conditions
    # on the electrostatic potential phi (note that periodic BCs on one but not the other does not make
    # sense and DECSKS will not accept such an input file; simulation will not start and InputErrors
    # will be broadcasted to the user)

    z.postpointmesh[k,:,:] = np.mod(z.postpointmesh[k,:,:], z.N)
    vz.postpointmesh[k,:,:] = vz.prepointmesh

    return f_k, f_old, Uf

def nonperiodic(f_k,
                f_old,
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

    NOTE: f_k is the distribution function contribution to all postpoints k.
    Only a symmetry boundary condition requires access to this object as it needs
    to make sure particles whose postpoints are exactly on the symmetry plane
    are factored in before replacing them in the reused f_old and Uf_old objects
    with their entering partners which are subsequently remapped as usual. This
    special consideration is described in github.com/dsirajud/../DECSKS-18 part 3 notebook.
    """
    # lower boundary
    f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['lower'])(
                                              f_old,
                                              Uf,
                                              z.postpointmesh[k,:,:],
                                              vz.prepointmesh,
                                              vz.prepointvaluemesh,
                                              z.N, vz.N, k, charge,
                                              sim_params,
                                              z, vz)

    # upper boundary
    f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['upper'])(
                                              f_old,
                                              Uf,
                                              z.postpointmesh[k,:,:],
                                              vz.prepointmesh,
                                              vz.prepointvaluemesh,
                                              z.N, vz.N, k, charge,
                                              sim_params,
                                              z, vz)


    # this orchestrator will be called if the distribution function boundary conditions are
    #
    #    (1) absorbing
    #    (2) collector
    #
    # these boundary conditions remove particles from the distribution function and eithier
    # absorbs them at walls with zero further influece (i.e. such a condition would be specified for
    # biased electrodes at a constant voltage), or their charge is collected and then asborbed.
    # the routine absorbs (zeroes out the density) and sets their postpoints to correspond
    # to such a wall where they are absorbed/collected. The remap procedure in the parent
    # function (convect_configuration or convect_velocity) then remaps these zeroed entries
    # to the gridpoint corresponding to the wall, which contribute zero. There is no numerical
    # need to change the velocity postpoints as the entries of the distribution themselves which
    # are zero do not affect the distribution function hereafter, so its velocity postpoints
    # are inconsequential.
    #
    vz.postpointmesh[k,:,:] = vz.prepointmesh

    return f_k, f_old, Uf

def symmetric(f_k,
              f_old,
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
    f_k, f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['lower'])(f_k,
                                                                                                   f_old,
                                                                                                   Uf,
                                                                                                   z.postpointmesh[k,:,:],
                                                                                                   vz.prepointmesh,
                                                                                                   vz.postpointmesh[k,:,:],
                                                                                                   vz.prepointvaluemesh,
                                                                                                   z.N, vz.N, k,
                                                                                                   charge,
                                                                                                   sim_params,
                                                                                                   z, vz)


    # upper boundary
    f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['upper'])(
                                              f_old,
                                              Uf,
                                              z.postpointmesh[k,:,:],
                                              vz.prepointmesh,
                                              vz.prepointvaluemesh,
                                              z.N, vz.N, k, charge,
                                              sim_params,
                                              z, vz)


    # vz.postpointmesh[k,:,:] to be set in boundaryconditions.symmetric_lower_boundary

    return f_k, f_old, Uf

# CYTHON METHODS
@cython.boundscheck(False)
def absorbing_lower_boundary(np.ndarray[DTYPE_t, ndim=2] f_old,
                              np.ndarray[DTYPE_t, ndim=2] Uf_old,
                              np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
                              np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
                              np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
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
                              np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
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
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
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

    #    j_inflection = sim_params['vx_inflection_gridpoint']
    if k == 0:
        for i in range(Nz):
            for j in range(Nvz):
                if zpostpointmesh[i,j] <= 0:
                    sigma_nk += -vzprepointvaluemesh[i,j] * (f_old[i,j] + Uf_old[i,j])
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = 0  # set postpoint at the absorber

    elif k == 1:
        for i in range(Nz):
            for j in range(Nvz):
                if zpostpointmesh[i,j] <= 0:
                    sigma_nk += vzprepointvaluemesh[i,j] * Uf_old[i,j]
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = 0  # set postpoint at the absorber


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
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
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

    #j_inflection = sim_params['vx_inflection_gridpoint']

    if k == 0:
        for i in range(Nz):
            for j in range(Nvz):
                if zpostpointmesh[i,j] >=  Nz - 1:
                    sigma_nk += vzprepointvaluemesh[i,j] * (f_old[i,j] - Uf_old[i,j])
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = Nz - 1  # set postpoint at the absorber

    elif k == 1:
        for i in range(Nz):
            for j in range(Nvz):
                if zpostpointmesh[i,j] >=  Nz - 1:
                    sigma_nk += vzprepointvaluemesh[i,j] * Uf_old[i,j]
                    f_old[i,j] = 0
                    Uf_old[i,j] = 0
                    zpostpointmesh[i,j] = Nz - 1  # set postpoint at the absorber


    sigma_nk *= charge
    sigma_nk *= vzwidth

    # update cumulative charge density
    sim_params['sigma'][z.str]['upper'] += sigma_nk
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old

@cython.boundscheck(False)
def symmetric_lower_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_k,
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzpostpointmesh,
        np.ndarray[DTYPE_t, ndim=2] vzprepointvaluemesh,
        int Nz, int Nvz, int k,
        int charge,
        sim_params,
        z, vz
        ):
    """

    This routine assumes a time substep is positive as the Vlasov-Poisson system are irreversible
    when boundaries are present (i.e. the upper boundary cannot be periodic/open in this setup)

    this boundary condition must not inherently must operate on a symmetric velocity grid,
    but it is more natural and consistent this way.

    For example, if an exiting particle at the symmetry plane has
    a velocity V1 < 0, it is the case that there exists an entering partner particle from the
    complement domain with oppositely directly velocity -V1 > 0. If a symmetric velocity
    grid, i.e. [avx, bvx] = [-V, V], is used then it is always possible to factor in both
    exiting (V1 < 0) density with entering partners (-V1 > 0) provided that V1 is inside [-V,V].

    If we have an asymmetric case where [avx, bvx] where avx != -bvx, then in the handling
    of the symmetry plane, there would be cases where such an exiting particle having V1 < 0 could
    its entering partner with velocity -V1 could not be represented (i.e. if avx < V1 is in the velocity grid,
    but bvx < -V1). There are two options for handling this, one of which is erroneous.

        (1) [Erroneous] The routine could try to muster some approximation to the scenario by
            taking the maximum value on the grid, but that would be unphysical.

        (2) [Correct] We could still argue that it would be appropriate means of handling such
            situations would be that if a particle with V1 < 0 exits, and if an entering particle
            would have -V1 > bvx, then to take the partner density to be zero and interpret the
            situation as the "window" of phase space we are looking at just does not show the entering partner.


    method (2) is permissible; however, we do not have current motivation to code this generality
    as the circumstances seem inferior to symmetric velocity grid setups. Instead, we currently
    advise users to set up symmetric velocity grids if a symmetry boundary condition is used.
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
