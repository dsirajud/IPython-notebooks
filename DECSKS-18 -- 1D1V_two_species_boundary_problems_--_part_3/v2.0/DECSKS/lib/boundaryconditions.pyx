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
             split_coeff,
             charge,
             k = 0
             ):
    """Applies periodic boundary conditions to postpointmesh on
    the active (variable being evolved) z.

    For such a scenario, the postpointmesh for the generalized
    velocity vz will remain the same as periodic boundary conditions
    will not produce any changes in vz.

    inputs:
    f_k -- (ndarray, ndim=2)
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

    return f_k, f_old, Uf

def nonperiodic(f_k,
                f_old,
                Uf,
                z,
                vz,
                sim_params,
                split_coeff,
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
    f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['lower'])(f_old, Uf,
                                                                                                  z.postpointmesh[k,:,:],
                                                                                                  vz.prepointmesh,
                                                                                                  z.N, vz.N, k, charge,
                                                                                                  sim_params,
                                                                                                  z, vz)

    # upper boundary
    f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['upper'])(f_old, Uf,
                                                                                                  z.postpointmesh[k,:,:],
                                                                                                  vz.prepointmesh,
                                                                                                  z.N, vz.N, k, charge,
                                                                                                  sim_params,
                                                                                                  z, vz)


    vz.postpointmesh[k,:,:] = vz.prepointmesh

    return f_k, f_old, Uf

def symmetric(f_k,
              f_old,
                Uf,
                z,
                vz,
                sim_params,
                split_coeff,
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
    # lower boundary, since the symmetric boundary condition orchestrator has been called, this means
    # the lower boundary condition routine will be lib.boundaryconditions.symmetric_lower_boundary
    f_k, f_old, Uf =  eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['lower'])(f_k,
                                                                                                         f_old, Uf,
                                                                                                         z.postpointmesh[k,:,:],
                                                                                                         vz.prepointmesh,
                                                                                                         z.N, vz.N, k,
                                                                                                         split_coeff,
                                                                                                         charge,
                                                                                                         sim_params,
                                                                                                         z, vz)

    # upper boundary
    f_old, Uf =  eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['upper'])(f_old, Uf,
                                                                                                    z.postpointmesh[k,:,:],
                                                                                                    vz.prepointmesh,
                                                                                                    z.N, vz.N, k, charge,
                                                                                                    sim_params,
                                                                                                    z, vz)



    # vz.postpointmesh[k,:,:] to be set in boundaryconditions.symmetric_lower_boundary
    # or boundaryconditions.symmtric_upper_boundary (if coded)

    return f_k, f_old, Uf

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
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        int Nz, int Nvz, int k, int charge,
        sim_params,
        z, vz
        ):
    # vars here are typed as C data types to minimize python interaction

    # for any such density packet whose postpoint predicts exceeding the lower boundary of
    # the domain, zero out (absorb) the density, and add contribution to total charge
    # density sigma at that boundary
    cdef int i, j
    cdef DTYPE_t sigma_n = 0 # at current time step n
    cdef DTYPE_t vzwidth = vz.width
    for i in range(Nz):
        for j in range(Nvz):
            if zpostpointmesh[i,j] <= 0:
                sigma_n += vzprepointmesh[i,j] * f_old[i,j]
                f_old[i,j] = 0
                Uf_old[i,j] = 0
                zpostpointmesh[i,j] = 0 # set postpoint at the absorber

    sigma_n *= charge
    sigma_n *= vzwidth

    # update cumulative charge density
    sim_params['sigma'][z.str]['lower'] += sigma_n
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old

@cython.boundscheck(False)
def collector_upper_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        int Nz, int Nvz, int k,
        int charge,
        sim_params,
        z, vz
        ):
    # vars here are typed as C data types to minimize python interaction

    # for any such density packet whose postpoint predicts exceeding the lower boundary of
    # the domain, zero out (absorb) the density, and add contribution to total charge
    # density sigma at that boundary
    cdef int i, j
    cdef DTYPE_t sigma_n = 0 # at current time step n
    cdef DTYPE_t vzwidth = vz.width
    for i in range(Nz):
        for j in range(Nvz):
            if zpostpointmesh[i,j] >=  Nz - 1:
                sigma_n += vzprepointmesh[i,j] * f_old[i,j]
                f_old[i,j] = 0
                Uf_old[i,j] = 0
                zpostpointmesh[i,j] = Nz - 1  # set postpoint at the absorber

    sigma_n *= charge
    sigma_n *= vzwidth

    # update cumulative charge density
    sim_params['sigma'][z.str]['upper'] += sigma_n
    z.postpointmesh[k,:,:] = zpostpointmesh

    return f_old, Uf_old

@cython.boundscheck(False)
def symmetric_lower_boundary(
        np.ndarray[DTYPE_t, ndim=2] f_k,
        np.ndarray[DTYPE_t, ndim=2] f_old,
        np.ndarray[DTYPE_t, ndim=2] Uf_old,
        np.ndarray[DTYPEINT_t, ndim=2] zpostpointmesh,
        np.ndarray[DTYPEINT_t, ndim=2] vzprepointmesh,
        int Nz, int Nvz, int k,
        float split_coeff,
        int charge,
        sim_params,
        z, vz
        ):

    cdef int i, j
    if k == 0:
        if Nvz % 2 == 0: # number of z gridpoints is even
            if split_coeff >= 0: # then vz < 0 particles can exit on left, 0 <= j <= Nvz // 2
                for j in range(0, Nvz / 2):
                    for i in range(Nz):
                        if z.postpointmesh[k,i,j] < 0:
                            # prepare for partner remap in lib.convect_configuration.remap_step
                            z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j]
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j]

                        elif z.postpointmesh[k,i,j] == 0:

                            # remap fraction of exiting density packet: f_k += f_old + Uf_old, U < 0
                            f_k += f_old[i,j]
                            f_k += Uf_old[i,j]

                            # prepare for partner remap in lib.convect_configuration.remap_step
                            # z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j] = 0 already
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j] # for subsequent remapping of the entering partner

            else: # split_coeff < 0: then vz > 0 particles can exit on left, Nvz // 2 <= j <= Nvz-1
                for j in range(Nvz / 2, Nvz):
                    for i in range(Nz):
                        if z.postpointmesh[k,i,j] < 0:
                            # prepare for partner remap in lib.convect_configuration.remap_step
                            z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j]
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j]

                        elif z.postpointmesh[k,i,j] == 0:

                            # remap fraction of exiting density packet: f_k += f_old + Uf_old, U < 0
                            f_k += f_old[i,j]
                            f_k += Uf_old[i,j]

                            # prepare for partner remap in lib.convect_configuration.remap_step
                            # z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j] = 0 already
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j] # for subsequent remapping of the entering partner

        else: # number of z gridpoints is odd (Nvz % 2 == 1), can skip over median index j; saves Nx loops
            if split_coeff >= 0: # then vz < 0 particles can exit on left, 0 <= j <= Nvz // 2
                for j in range(0, Nvz / 2):
                    for i in range(Nz):
                        if z.postpointmesh[k,i,j] < 0:
                            # prepare for partner remap in lib.convect_configuration.remap_step
                            z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j]
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j]

                        elif z.postpointmesh[k,i,j] == 0:

                            # remap fraction of exiting density packet: f_k += f_old + Uf_old, U < 0
                            f_k += f_old[i,j]
                            f_k += Uf_old[i,j]

                            # prepare for partner remap in lib.convect_configuration.remap_step
                            # z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j] = 0 already
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j] # for subsequent remapping of the entering partner

            else: # split_coeff < 0: then vz > 0 particles can exit on left, Nvz // 2 <= j <= Nvz-1
                for j in range(Nvz / 2 + 1, Nvz):
                    for i in range(Nz):
                        if z.postpointmesh[k,i,j] < 0:
                            # prepare for partner remap in lib.convect_configuration.remap_step
                            z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j]
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j]

                        elif z.postpointmesh[k,i,j] == 0:

                            # remap fraction of exiting density packet: f_k += f_old + Uf_old, U < 0
                            f_k += f_old[i,j]
                            f_k += Uf_old[i,j]

                            # prepare for partner remap in lib.convect_configuration.remap_step
                            # z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j] = 0 already
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j] # for subsequent remapping of the entering partner


    if k == 1:
        if Nvz % 2 == 0: # number of z gridpoints is even
            if split_coeff >= 0: # then vz < 0 particles can exit on left, 0 <= j <= Nvz // 2
                for j in range(0, Nvz / 2):
                    for i in range(Nz):
                        if z.postpointmesh[k,i,j] < 0:
                            # prepare for partner remap in lib.convect_configuration.remap_step
                            z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j]
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vz.postpointmesh[k,i,j]
                            Uf_old[i,j] = -Uf_old[i,j]

                        elif z.postpointmesh[k,i,j] == 0:

                            # remap fraction of exiting density packet: f_k += f_old + Uf_old, U < 0
                            f_k -= Uf_old[i,j]

                            # prepare for partner remap in lib.convect_configuration.remap_step
                            # z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j] = 0 already
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j] # for subsequent remapping of the entering partner

            else: # split_coeff < 0: then vz > 0 particles can exit on left, Nvz // 2 <= j <= Nvz-1
                for j in range(Nvz / 2, Nvz):
                    for i in range(Nz):
                        if z.postpointmesh[k,i,j] < 0:
                            # prepare for partner remap in lib.convect_configuration.remap_step
                            z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j]
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vz.postpointmesh[k,i,j]
                            Uf_old[i,j] = -Uf_old[i,j]

                        elif z.postpointmesh[k,i,j] == 0:

                            # remap fraction of exiting density packet: f_k += f_old + Uf_old, U < 0
                            f_k -= Uf_old[i,j]

                            # prepare for partner remap in lib.convect_configuration.remap_step
                            # z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j] = 0 already
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j] # for subsequent remapping of the entering partner

        else: # number of z gridpoints is odd (Nvz % 2 == 1), can skip over median index j; saves Nx loops
            if split_coeff >= 0: # then vz < 0 particles can exit on left, 0 <= j <= Nvz // 2
                for j in range(0, Nvz / 2):
                    for i in range(Nz):
                        if z.postpointmesh[k,i,j] < 0:
                            # prepare for partner remap in lib.convect_configuration.remap_step
                            z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j]
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j]

                        elif z.postpointmesh[k,i,j] == 0:

                            # remap fraction of exiting density packet: f_k += f_old + Uf_old, U < 0
                            f_k -= Uf_old[i,j]

                            # prepare for partner remap in lib.convect_configuration.remap_step
                            # z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j] = 0 already
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vz.postpointmesh[k,i,j]
                            Uf_old[i,j] = -Uf_old[i,j] # for subsequent remapping of the entering partner

            else: # split_coeff < 0: then vz > 0 particles can exit on left, Nvz // 2 <= j <= Nvz-1
                for j in range(Nvz / 2 + 1, Nvz):
                    for i in range(Nz):
                        if z.postpointmesh[k,i,j] < 0:
                            # prepare for partner remap in lib.convect_configuration.remap_step
                            z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j]
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j]

                        elif z.postpointmesh[k,i,j] == 0:

                            # remap fraction of exiting density packet: f_k += f_old + Uf_old, U < 0
                            f_k -= Uf_old[i,j]

                            # prepare for partner remap in lib.convect_configuration.remap_step
                            # z.postpointmesh[k,i,j] = -z.postpointmesh[k,i,j] = 0 already
                            vz.postpointmesh[k,i,j] = Nvz - 1 - vzprepointmesh[i,j]
                            Uf_old[i,j] = -Uf_old[i,j] # for subsequent remapping of the entering partner


    return f_k, f_old, Uf_old
