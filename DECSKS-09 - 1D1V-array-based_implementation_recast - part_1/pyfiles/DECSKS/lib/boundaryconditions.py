import numpy as np

def periodic(sim_params, z, Uf = None):
    """Applies periodic boundary conditions to
    postpoints

    inputs:
    z -- (instance) phase space variable. Used to
         update attribute z.postpointmesh
         (2D mesh of advected gridpoints, not values)

    outputs:
    z.postpoints -- (attr) update to z.postpoints
                    attribute per periodic BCs
    """

    # there isn't a reason for these if statements
    # should just do the latter, and then for the ballistic
    # push, use k[0,:,:]

    # if calling for Lagrangian integral push
    if Uf is None:
        z.postpointmesh = np.mod(z.postpointmesh, z.N)
        return z.postpointmesh

    # if calling for remap indices
    if Uf is not None:
        # remap index meshes, k1[i,j], k2[i,j] correspond to
        # mapped to indices of prepointmesh[i,j]

        kshape = [2]
        for dim in sim_params['dims']:
            kshape.append(dim)

        k = np.zeros(kshape)
        k[0,:,:] = z.postpointmesh
        k[1,:,:] = np.where(Uf >= 0, np.mod(z.postpointmesh + 1 , z.N),
                      np.mod(z.postpointmesh - 1, z.N))

        return k

def periodic_old(z, i = None, Uf = None):
    """Applies periodic boundary conditions to
    postpoints

    inputs:
    z -- (instance) phase space variable
    i -- (int)
    Uf -- (ndarray, ndim=1) normalized fluxes
          for z.MCs

    outputs:
    k1, k2 -- (int) postpoint indices per
              periodic BCs
    """

    if i is None and Uf is None:
        # whole vector z.postpoints passed
        z.postpoints = np.mod(z.postpoints, z.N)
        return z.postpoints

    else:
        if Uf[i] > 0:
            k1 = z.postpoints[i]
            k2 = np.mod(z.postpoints[i] + 1, z.N)
        elif Uf[i] < 0:
            k1 = z.postpoints[i]
            k2 = np.mod(z.postpoints[i] - 1, z.N)
        else:
            # dummy entries, will store in summary
            # array izf as float('nan')
            k1 = None
            k2 = None
        return k1, k2
