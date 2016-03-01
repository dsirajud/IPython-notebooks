import numpy as np

def periodic(f_old,
             Uf,
             z,
             sim_params, # used in boundaryconditions.nonperiodic
             k = None
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
    z.postpointmesh = np.mod(z.postpointmesh, z.N)

    return f_old, Uf, z

def periodize_postpointmesh(zpostpointmesh, zN):
    """Applies periodic boundary conditions to
    postpointmesh[k,:,:]

    inputs:
    z.postpointmesh -- (ndarray, ndim=2), shape = (x.N, vx.N)

    outputs:
    z.postpointmesh -- (ndarray, ndim=2), shape = (x.N, vx.N)
                        periodic BCs applied

    """
    zpostpointmesh = np.mod(zpostpointmesh, zN)

    return zpostpointmesh

def nonperiodic(f_old,
                Uf,
                z,
                sim_params,
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
                           '_lower_boundary')(f_old, Uf,
                                              z.postpointmesh[k,:,:])

    # upper boundary
    f_old, Uf = eval(sim_params['BC'][z.str]['upper'] +
                           '_upper_boundary')(f_old, Uf,
                                              z.postpointmesh[k,:,:], z)

    # since the relevant entries of f_old and Uf that exit the domain
    # are zeroed out, in order to have a clean addition as before
    # we map their postpoints to their corresponding periodic BC locations
    # so that there are no two shared postpoints by construction

    z.postpointmesh[k,:,:] = periodize_postpointmesh(z.postpointmesh[k,:,:], z.N)
    # note that if there were shared postpoints, then lib.convect.remap_assignment
    # would not allocate the correct densities since it is a matrix sum
    # rather than a (slow) loop (where we could have used +=). For example
    # should more than one prepoint share a common postpointss, only one
    # cell's density would be allocated to the postpoint, the rest would be
    # overwritten, not incrementally summed

    # note the other function, periodic, periodizes z.postpointmesh.shape = (2, z.N, vz.N)
    # here, we only want to periodize the postpointmesh pertaining to the index 'nearest'
    # or 'contiguous' so that we do not tarnish the postpointmesh when passing through 'nearest'
    # so that 'contiguous' is unphysically periodized, hence most boundaries are evaded.
    return f_old, Uf, z

def absorbing_lower_boundary(f_old, Uf, zpostpointmesh):

    f_old = np.where(zpostpointmesh <= 0, 0, f_old)
    Uf = np.where(zpostpointmesh <= 0, 0, Uf)

    return f_old, Uf

def absorbing_upper_boundary(f_old, Uf, zpostpointmesh, z):

    f_old = np.where(zpostpointmesh >= z.N, 0, f_old)
    Uf = np.where(zpostpointmesh >= z.N, 0, Uf)

    return f_old, Uf
