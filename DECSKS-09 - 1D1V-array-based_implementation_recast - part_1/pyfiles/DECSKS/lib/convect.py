import DECSKS
import numpy as np
import numpy.ma as ma

def scheme(
    f_initial,
    n,
    sim_params,
    z,
    vz
    ):
    """Solve a 1D advection (in z) equation by Convected Scheme

    inputs:
    f_old -- (ndarray, ndim = 1 or 2) f[n-1,:] or f[n-1,:,:] if first
             substep in a splitting algorithm or if none
             else, f[n,:] or f[n,:]
    z -- (instance) phase space variable
    n -- (int) time step
    sim_params -- (dict) simulation parameters

    outputs:
    f_final -- (ndarray, ndim = 1 or 2) f[n,:] or f[n,:,:] updated
               after all steps have been completed
    """

    # (0) INITIALIZE FINAL DENSITY CONTAINER AND EXTRACT EVOLVED GRID
    f_final = np.zeros_like(f_initial)
    f_initial = DECSKS.lib.convect.extract_active_grid(sim_params, f_initial)

    # (1) PUSH MOVING CELLS an integer number of cells
    z, k = advection_step(f_initial, z, vz)

    # (*) APPLY BOUNDARY CONDITIONS
    z.postpointmesh = DECSKS.lib.boundaryconditions.periodic(z.postpointmesh, z.N)
    k = DECSKS.lib.boundaryconditions.periodic(k, z.N)

    # (2) REMAP DENSITY TO GRID
    f_remapped = remap_step(
                       sim_params,
                       f_initial,
                       k,
                       n,
                       z,
                       vz
                       )

    # (3) COLLISION STEP (NOT YET IMPLEMENTED)
    # f_new = DECSKS.lib.collisions.collisiontype(f_old, z, n)

    # (4) RETURN FINAL DESTINY (density*)
    f_final = finalize_density(sim_params, f_remapped, f_final)

    return f_final
#---------------------------------------------------------------------------  #
def advection_step(f_old, z, vz):
    """Pushes each z.prepointmesh (index) value by the advection *along each
    column* (i.e. constant vz.prepointvaluemesh)
    as prescribed by its generalized velocity vz.prepointmeshvalues according to

        vz.prepointmeshvalues * dt / z.width = CFL.numbers

    the integer part, CFL.int, is pushed here. The residual fraction
    CFL.frac is remapped according to remap_step(*args) in step (2) in
    the orchestrator routine, scheme(*args), above.

    inputs:
    sim_params -- (dict) simulation parameters
    z -- (instance) phase space variable equipped with (among other attributes)

            z.prepointmesh -- (ndarray, ndim=2) initial [i,j] prepoint indices of all MCs
            z.MCs -- (ndarray, ndim = 2) final (non-integral) [i'',j'']
                     postpoint indices for all MCs after push

    outputs:
    z -- (instance) phase space variable updated with postpointmesh attribute (BCs not applied)
    k -- (ndarray, ndim=2) matrix containing each pair of postpoints for each prepoint [i,j]

    NOTE: boundary conditions not applied at this point to z, k
    """
    z.postpointmesh = z.prepointmesh + z.CFL.int
    k = generate_postpoint_mesh(f_old, z, vz)

    return z, k

def generate_postpoint_mesh(
                            f_old,
                            z,
                            vz):
    """
    generates the postpoint mapping matrix k --
                k[0,:,:] gives the nearest neighbor gridpoint for
                         each MC's postposition that originated at [i,j]

                k[1,:,:] gives the contiguous gridpoint for
                         each MC's postposition that originated at [i,j]
                         (depends on the sign of vz)

                k.shape = (2, f_old.shape[0], f_old.shape[1])
                          which is referenced from sim_params['dims']

                NOTE: k[r,i,j] for r = 0 or 1 gives the postpoint
                      of the advecting index, it does not give
                      the relevant postpoint duple in phase space:

                          (i + k[r,i,j], j) [for x advection]

                      or

                          (i, j + k[r,i,j]) [for v advection]

                      but instead just gives that value k[r,i,j].
                      It then makes sense to write statements such as

                          k[0,i,j] + 1 is a contiguous gridpoint to
                                       the index k[0,i,j] in the positive
                                       direction

                          k[0,i,j] - 1 is a contiguous gridpoint to
                                       the index k[0,i,j] in the negative
                                       direction

                      without risk of inadvertently updating the opposite
                      (non-advecting) gridpoint

    inputs:
    sim_params -- (dict) simulation parameters


    outputs:


    """
    k = np.zeros([2, f_old.shape[0], f_old.shape[1]])
    k[0,:,:] = z.postpointmesh

    # if velocity > 0, MC remapped to k[0,i,j] and (k[0,i,j] + 1)
    # elif velocity < 0, MC is remapped to k[0,i,j] and (k[0,i,j] - 1)
    k[1,:,:] = np.where(vz.prepointvaluemesh >= 0, z.postpointmesh + 1,
                  z.postpointmesh - 1)

    k = np.array(k, dtype = int)

    return k

def remap_step(
        sim_params,
        f_old,
        k,
        n,
        z,
        vz
        ):
    """Orchestrates remapping of all advected moving cells to the grid,
    the remapping rule is in the called method remap_assignment(*args)

    inputs:
    sim_params -- (dict) simulation parameters
    f_old -- (ndarray, dim=2) density from previous time step, full grid
    z.CFL -- (instance) contains z.CFL.frac attribute = (ndarray, ndim=2)
           which gives the fractional CFL number for every MC
           originating at prepoint duples {(i,j)} --> {(i + z.CFL.frac, j)} or
           {(i,j)} --> {(i, j + z.CFL.frac)} depending on which variable is being
           advected,

           i.e. these values are the distance normalized by the cell width to the
           nearest index (z.CFL.int) defined with the convention that z.CFL.int
           is rounded towards index zero. (Note, the user should not be tempted to
           use the Python floor function as it rounds towards -infinity, not zero.
           This is a perpetual subject of controversy in Python.
           The issue has been long resolved by the BDFL, the definition of floor
           (same as //, also same as / in python 3.1+) exists to preserve the simple
           arithmetic reconstruction of a divided number using /, %).

    n  -- (int) current time step
    z -- (instance) phase space variable
    vz -- (instance) generalized velocity variable for phase space variable z

    outputs:
    f_new -- (ndarray, dim=2) density with all MCs remapped at final postpoints
                     according to remap rule

    NOTE: we carry assemble the new density as a sum of two
    2D matrices f_k1 and f_k2 which themselves are the associated
    mappings to the postpoints k1 (nearest neighbor grid index to
    exact non-integral postpoint), and k2 (contiguous grid index
    to k1 based on the sign of advection)

    The reason a partition is necessary is that a mapping such as

        f_new[k, vz.prepointmesh] += f_old

    does not increment each element, but instead overwrites leaving
    only the final assignment to a particular gridpoint rather than
    adds on top of it. Thus, the separation of k1 and k2 into
    different containers f_k1 and f_k2 permit non-overlapping
    assignments in each individually, and their sum gives the
    appropriate total density f_new.
    """

    # compute high order fluxes
    Uf = flux(
        sim_params,
        f_old,
        z, vz
        )

    # remap to nearest neighbor cell center
    f_k1 = remap_assignment(
        f_old,
        Uf,
        k[0,:,:],
        z,
        vz,
        index = 'nearest'    # remaps to k1
        )

    # remap to contiguous cell center
    f_k2 = remap_assignment(
        f_old,
        Uf,
        k[1,:,:],
        z,
        vz,
        index = 'contiguous'    # remaps to k2
        )

    f_new = f_k1 + f_k2

    # global check on density conservation
    DECSKS.lib.density.global_conservation_check(sim_params, f_new, n)

    return f_new
# ........................................................................... #
def flux(
        sim_params,
        f_old,
        z, vz
        ):
    """Computes fluxes Uf for all z.prepointmesh

    inputs:
    f_old -- (ndarray, dim=2) density from previous time (sub)step
    CFL -- (instance) CFL number
    z -- (instance) phase space variable
    sim_params -- (dict) simulation parameters

    outputs:
    Uf -- (ndarray, dim=2) Normalized fluxes for every z.prepoints[i]

           where for each [i,j]:

               Uf[i,j] = sum c[q,j]*d[q,i,j] over q = 0, 1, ... N-1
    """

    c = DECSKS.lib.HOC.correctors(sim_params, z, vz)

    # container for derivative coefficients, d[q,i,j], for q-th deriv at [i,j]
    d = np.zeros([sim_params['N'], z.N, vz.N])

    # zeroeth derivative is f_old itself
    d[0,:,:] = f_old

    # evaluate derivatives q = 1, 2, ... N-1

    # calls lib.derivatives.fd or lib.derivatives.fourier based on the
    # HOC specified in etc/params.dat, sim_params['derivative_method']
    # contains the related function handle as a string
    d = eval(sim_params['derivative_method'])(d, f_old, z, sim_params)

    # compute high order fluxes column-by-column
    Uf = np.zeros(f_old.shape)
    for j in range(vz.N):
        Uf[:,j] = c[:,j].dot(d[:,:,j])

    # enforce flux limiter to ensure positivity and restrict numerical overflow
    Uf = flux_limiter(f_old, Uf, z)

    return Uf

def flux_limiter(
        f_old,
        Uf,
        z):
    """Applies the flux (numerical and positivity) limiter to the
    ith nominal flux matrix Uf = (z1.N, z2.N)

    inputs:
    f_old -- (ndarray, dim=2) density from previous time step
    CFL -- (instance) Courant numbers characterizing MC convection
    Uf -- (ndarray, dim=2) the normalized flux (high order or not) used in the CS update

    outputs:
    Uf -- (ndarray, dim=2) final normalized flux for MC originating at prepoints[i,j]
            after numerical and positivity limiter has been applied


    Note: for a 768 x 1536 grid

          the masked implementation here takes about 285 ms
          the looped implementation over all i,j takes 3.85 s for the whole grid

          for a 1536 x 3072 grid

          the masked implementation here takes about 1.44 s
          the looped implementation over all i,j takes 17.9 s for the whole grid

          i.e. the computational savings is at least a factor of 10
    """
    # local masked array (ma) copy of CFL.frac used to leave CFL.frac unchanged
    Uf_ma = ma.array(z.CFL.frac)

    # mask negative values, keep positives
    Uf_ma[z.CFL.frac < 0] = ma.masked

    # assign the mask to a zero matrix of same size for pairwise ma.minimum/maximum operations below
    zeros = ma.zeros(Uf.shape)
    zeros.mask = Uf_ma.mask

    # operating only on positive values (negatives masked out)
    Uf_pos = ma.minimum(ma.maximum(zeros, Uf), 1.0*f_old)

    # mask positives, keep negatives
    Uf_ma.mask = np.logical_not(Uf_ma.mask)
    zeros.mask = Uf_ma.mask

    #operating only on negative values
    Uf_neg = ma.minimum(ma.maximum(-1.0*f_old, Uf), zeros)

    # combine masked values in a single matrix Uf_final
    Uf = np.where(Uf_neg.mask == False, Uf_neg.data, Uf_pos.data)

    return Uf

def remap_assignment(
        f_old,
        Uf,
        k,
        z,
        vz,
        index = 'nearest'
        ):
    """Remaps the z.N MCs to Eulerian mesh at indices k1 and k2

    inputs:
    f_old -- (ndarray, dim=1) density from previous time step
    Uf -- (ndarray, dim=1) the normalized flux used in the CS update
    z -- (instance) phase space variable whose MCs are convected
    n -- (int) current time step

    outputs:
    f_new -- (ndarray, dim=1) density with f_old remapped to
             according to the mapping

            [z.prepointmesh, vz.prepointmesh] --> [k, vz.prepointmesh]

    """
    mask_neg =  (z.CFL.frac < 0)
    mask_pos = np.logical_not(mask_neg)

    f_old_ma = ma.array(f_old)
    f_pos, f_neg = ma.zeros(f_old.shape), ma.zeros(f_old.shape)

    Uf_ma = ma.array(Uf)

    if index == 'nearest':
        # mask out negative values
        f_old_ma.mask = mask_neg
        Uf_ma.mask = mask_neg

        f_pos[ k, vz.prepointmesh ] = f_old_ma - Uf_ma

        # mask out all positive values
        f_old_ma.mask = mask_pos
        Uf_ma.mask = mask_pos

        f_neg[ k, vz.prepointmesh ] = f_old_ma + Uf_ma

    elif index == 'contiguous':
        # mask out negative values
        f_old_ma.mask = mask_neg
        Uf_ma.mask = mask_neg

        f_pos[ k, vz.prepointmesh ] = Uf_ma

        # mask out all positive values
        f_old_ma.mask = mask_pos
        Uf_ma.mask = mask_pos

        f_neg[ k, vz.prepointmesh ] = -Uf_ma

    # consolidate data, "for every entry where mask_neg is False [positi]"
    f_new = np.where(mask_neg == True, f_neg.data, f_pos.data)

    return f_new

def extract_active_grid(sim_params, f_total_grid):
    """We evolve the density from the previous time step, f_old
    only on the gridpoints that are 'active' (cf. DECSKS-09 notebook)
    We distinguish, e.g. in 1D, the two attributes of a phase space
    variable z:

        z.Ngridpoints -- (int) total number of gridpoints
        z.N           -- (int) total number of 'active' gridpoints

    The total grid indices  : 0, 1, ... , z.Ngridpoints - 1
    The active grid indices : 0, 1, ... , z.N - 1

    For all but periodic boundary conditions (PBCs), these are the same.
    That is, for periodic boundary conditions (PBCs):

        z.N = z.Ngridpoints - 1

    so we evolve f_old[:z.N] -> f_new[:z.N]

    and we complete by the density by periodicity:

        f_new[z.Ngridpoints - 1] = f_new[0]

    for all other BCs: z.N = z.Ngridpoints and this function has no
    effect.

    inputs:
    f_total_grid -- (ndarray, ndim=2) 2D density constituting total grid

    outputs:
    f_active_grid -- (ndarray, ndim=2) 2D density containing all
                     active gridpoints

    """

    f_active_grid = f_total_grid[0:sim_params['active_dims'][0], 0:sim_params['active_dims'][1]]

    return f_active_grid

def finalize_density(sim_params, f_remapped, f_final):
    """
    returns a final density. For all but PBCs, f_new = f_final since
    our setup is such that (z1.N, z2.N) moving cells are evolved.

    The bookkeeping is such that, e.g. in 1D

          non-periodic BCs : z.N = z.Ngridpoints

          periodic BCs     : z.N = z.Ngridpoints - 1
                             f[z.Ngridpoints-1] = f[0] by periodicity

    We use the associated generalization to two dimensions, e.g.

          non-periodic BCs: z1.N = z1.Ngridpoints
                            z2.N = z2.Ngridpoints

                            hence f_final = f_new

          periodic BCs     : z.N = z.Ngridpoints - 1
                             f_final[z1.Ngridpoints - 1, :] = f_new[0,:]
                             f_final[:, z2.Ngridpoints - 1] = f_new[:,0]
  """
    # TODO currently assuming that both dimensions are periodic, need to generalize
    # this in the future

    # assign all active grid points to grid values on f_final
    f_final[:f_remapped.shape[0], :f_remapped.shape[1]] = f_remapped

    f_final[f_remapped.shape[0]+1,:] = f_remapped[0,:]
    f_final[:, f_remapped.shape[1] + 1] = f_remapped[:,0]


    return f_final
