import numpy as np
import DECSKS

def periodic(f_old,
             Uf,
             z,
             vz,
             sim_params, # used in boundaryconditions.nonperiodic
             charge,
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
    f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['lower'])(f_old, Uf,
                                                                                                   z.postpointmesh[k,:,:], k,
                                                                                                   z, vz,
                                                                                                   sim_params, charge)


    f_old, Uf = eval(sim_params['distribution_function_boundarycondition_handle'][z.str]['upper'])(f_old, Uf,
                                                                                                   z.postpointmesh[k,:,:], k,
                                                                                                   z, vz,
                                                                                                   sim_params, charge)



    # since the relevant entries of f_old and Uf that exit the domain
    # are zeroed out, in order to have a clean addition as before
    # we map their postpoints to their corresponding periodic BC locations
    # so that there are no two shared postpoints by construction

    # map z.postpointmesh to periodic locations
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

def absorbing_lower_boundary(f_old, Uf, zpostpointmesh, k, z, vz, sim_params, charge):

    f_old = np.where(zpostpointmesh <= 0, 0, f_old)
    Uf = np.where(zpostpointmesh <= 0, 0, Uf)

    return f_old, Uf

def absorbing_upper_boundary(f_old, Uf, zpostpointmesh, k, z, vz, sim_params, charge):

    f_old = np.where(zpostpointmesh >= z.N-1, 0, f_old)
    Uf = np.where(zpostpointmesh >= z.N-1, 0, Uf)

    return f_old, Uf

def collector_lower_boundary(f_old, Uf, zpostpointmesh, k, z, vz, sim_params, charge):

    # this discriminates vx vs. x, as the boundary condition function handle
    # sim_params['BC'][z.str]['lower' or 'upper'] for z.str = 'vx' is never set to 'charge_collection'
    # in lib.read

    # collect density packets that reach the boundary (or beyond)

    Uf_absorbed = np.where(zpostpointmesh <= 0, Uf, 0)
    if k == 0: # if nearest grid point
        f_absorbed = np.where(zpostpointmesh <= 0, f_old, 0)

        # the fraction (1 + U[i,j]) of f_absorbed[i,j] is deposited at the wall (written below as (f + Uf))
        # any vz.prepointvaluemesh[i,j] which pushed a particle to this edge is negative
        # any Uf_absorbed which corresponds to the normalized flux of particles at the edge is also negative
        # the minus sign in front of vz.prepointvaluemesh acts as the absolute value operator
        sigma_nk = charge * np.sum(-vz.prepointvaluemesh * (f_absorbed + Uf_absorbed) * vz.width)

    elif k == 1:

        # the fraction (-U[i,j]) of f_absorbed[i,j] is deposited at the wall (written below as (-Uf))

        # minus signs have been cancelled here, the raw version of this might be thought of as
        # sigma_nk = charge * np.sum((-vz.prepointvaluemesh) * (-Uf_absorbed) * vz.width)
        # where Uf_absorbed < 0 as well.

        sigma_nk = charge * np.sum(vz.prepointvaluemesh * (Uf_absorbed) * vz.width)

    # update global variable, sigma
    sim_params['sigma'][z.str]['lower'] += sigma_nk
    # remove the exiting particles from f, Uf to prep for subsequent remapping (will remap with zero contribution)
    f_old, Uf = absorbing_lower_boundary(f_old, Uf, zpostpointmesh, k, z, vz, sim_params, charge)

    return f_old, Uf

def collector_upper_boundary(f_old, Uf, zpostpointmesh, k, z, vz, sim_params, charge):

    # this discriminates vx vs. x, as the boundary condition function handle
    # sim_params['BC'][z.str]['lower' or 'upper'] for z.str = 'vx' is never set to 'charge_collection'
    # in lib.read

    # collect density packets that reach the boundary (or beyond)

    Uf_absorbed = np.where(zpostpointmesh >= z.N-1, Uf, 0)
    if k == 0: # if nearest grid point
        f_absorbed = np.where(zpostpointmesh >= z.N-1, f_old, 0)

        # the fraction (1 - U[i,j]) of f_absorbed[i,j] is deposited at the wall (written below as (f - Uf))
        # any vz.prepointvaluemesh[i,j] which pushed a particle to this edge is positive
        # any Uf_absorbed which corresponds to the normalized flux of particles at the edge is also positive
        sigma_nk = charge * np.sum(vz.prepointvaluemesh * (f_absorbed - Uf_absorbed) * vz.width)

    elif k == 1:

        # the fraction (U[i,j]) of f_absorbed[i,j] is deposited at the wall (written below as (Uf))

        sigma_nk = charge * np.sum(vz.prepointvaluemesh * (Uf_absorbed) * vz.width)

    # update global variable, sigma
    sim_params['sigma'][z.str]['upper'] += sigma_nk
    # remove the exiting particles from f, Uf to prep for subsequent remapping (will remap with zero contribution)
    f_old, Uf = absorbing_upper_boundary(f_old, Uf, zpostpointmesh, k, z, vz, sim_params, charge)

    return f_old, Uf
