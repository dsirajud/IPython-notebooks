import numpy as np

# Python methods

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
                           '_lower_boundary')(f_old, Uf,
                                              z.postpointmesh[k,:,:], z, vz,
                                              sim_params, charge)

    # upper boundary
    f_old, Uf = eval(sim_params['BC'][z.str]['upper'] +
                           '_upper_boundary')(f_old, Uf,
                                              z.postpointmesh[k,:,:], z, vz,
                                              sim_params, charge)

    z.postpointmesh[k,:,:] = np.mod(z.postpointmesh[k,:,:], z.N)

    return f_old, Uf


# Cython methods

def symmetric_lower_boundary(f_old, Uf, zpostpointmesh, z, vz, sim_params, charge):
    f_entering = np.where(zpostpointmesh < 0, f_old, 0) # = f_exiting
    Uf_entering = np.where(zpostpointmesh < 0, -Uf, Uf)

def absorbing_lower_boundary(f_old, Uf, zpostpointmesh, z, vz, sim_params, charge):

    f_old = np.where(zpostpointmesh <= 0, 0, f_old)
    Uf = np.where(zpostpointmesh <= 0, 0, Uf)

    return f_old, Uf

def absorbing_upper_boundary(f_old, Uf, zpostpointmesh, z, vz, sim_params, charge):

    f_old = np.where(zpostpointmesh >= z.N, 0, f_old)
    Uf = np.where(zpostpointmesh >= z.N, 0, Uf)

    return f_old, Uf

def charge_collection_lower_boundary(f_old, Uf, zpostpointmesh, z, vz, sim_params, charge):

    # this discriminates vx vs. x, as the boundary condition function handle
    # sim_params['BC'][z.str]['lower' or 'upper'] for z.str = 'vx' is never set to 'charge_collection'
    # in lib.read
    f_absorbed = np.where(z.postpointmesh <= 0, f_old, 0)
    sigma_n = np.sum(vz.prepointmesh * f_absorbed * vz.width)

    # passed by reference, original value is modified, no need for explicit return
    sim_params['sigma'][z.str]['lower'] = \
      sim_params['sigma'][z.str]['lower'] + charge*sigma_n

    f_old, Uf = absorbing_lower_boundary(f_old, Uf, zpostpointmesh, z, vz, sim_params, charge)

    return f_old, Uf

def charge_collection_upper_boundary(f_old, Uf, zpostpointmesh, z, vz, sim_params, charge):

    # this discriminates vx vs. x, as the boundary condition function handle
    # sim_params['BC'][z.str]['lower' or 'upper'] for z.str = 'vx' is never set to 'charge_collection'
    # in lib.read
    f_absorbed = np.where(z.postpointmesh <= 0, f_old, 0)
    sigma_n = np.sum(vz.prepointmesh * f_absorbed * vz.width)

    # passed by reference, no need for explicit return
    sim_params['sigma'][z.str]['upper'] = \
      sim_params['sigma'][z.str]['upper'] + charge*sigma_n

    f_old, Uf = absorbing_upper_boundary(f_old, Uf, zpostpointmesh, z, vz, sim_params, charge)

    return f_old, Uf
