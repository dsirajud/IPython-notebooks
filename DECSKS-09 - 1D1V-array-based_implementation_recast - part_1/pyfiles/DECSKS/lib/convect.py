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
    """Solves a collection of 1D advection (in z) equations by convected scheme
    and stacks the results in a 2D matrix

    inputs:
    f_initial -- (ndarray, ndim = 2) f[n-1,:,:] if first
             substep in a splitting algorithm or if none
             else, f[n,:,:]
    n -- (int) time step
    sim_params -- (dict) simulation parameters
    z -- (instance) phase space variable
    vz -- (instance) generalized velocity for z


    outputs:
    f_final -- (ndarray, ndim = 1 or 2) f[n,:] or f[n,:,:] updated
               after all steps have been completed
    """

    # (0) INITIALIZE FINAL DENSITY CONTAINER AND EXTRACT EVOLVED GRID
    f_final = np.zeros(f_initial.shape)
    f_initial = DECSKS.lib.domain.extract_active_grid(f_initial, z, sim_params)

    if z.str[0] == 'v':
        f_initial, z, vz = DECSKS.lib.domain.velocity_advection_prep(f_final, f_initial, z, vz)

    # (1) PUSH MOVING CELLS an integer number of cells
    z.postpointmesh = advection_step(z)

    # (*) APPLY BOUNDARY CONDITIONS
    z.postpointmesh = DECSKS.lib.boundaryconditions.periodic(z.postpointmesh, z.N)

    # (2) REMAP DENSITY TO GRID
    f_remapped = remap_step(
                       sim_params,
                       f_initial,
                       n,
                       z,
                       vz
                       )

    # (3) COLLISION STEP (NOT YET IMPLEMENTED)
    # f_new = DECSKS.lib.collisions.collisiontype(f_old, z, n)

    # (4) RETURN FINAL DESTINY (density*)
    f_final = finalize_density(sim_params, f_remapped, f_final, z, vz)

    return f_final
#---------------------------------------------------------------------------  #
def advection_step(z):
    """Pushes each z.prepointmesh (index) value by the advection *along each
    column j* (i.e. constant vz.prepointvaluemesh)
    as prescribed by its generalized velocity vz.prepointmeshvalues[:,j].

    This is computed in one calculation by:

        vz.prepointmeshvalues * dt / z.width = CFL.numbers

    the integer part, CFL.int, is pushed here. The residual fraction
    CFL.frac is remapped according to remap_step(*args) in step (2) in
    the orchestrator routine, scheme(*args), above.

    inputs:
    z -- (instance) phase space variable equipped with (among other attributes)

            z.prepointmesh -- (ndarray, ndim=2) initial [i,j] prepoint indices of all MCs

            z.CFL -- (instance) contains CFL numbers, int and frac parts
                z.CFL.numbers -- (ndarray, ndim=2, dtype=float64)
                z.CFL.int -- (ndarray, ndim=2, dtype=int) integer part of CFL numbers
                z.CFL.frac -- (ndarray, ndim=2, dtype=float64) fractional part

    outputs:
    z -- (instance) phase space variable updated with postpointmesh attribute (BCs not applied)

        updated attr:
        z.postpointmesh -- (ndarray, ndim=3, dtype=int), shape = (2, x.N, v.N) matrix containing each pair
                           k[:,i,j] of postpoints for each prepoint [i,j]


    NOTE: Boundary conditions NOT applied at this point

    USAGE NOTE: the map z.postpointmesh[k,i,j] for k = 0 or 1 gives the postpoint
      of the advecting index, it does not give the relevant postpoint duple
      in phase space:

          (z.postpointmesh[k,i,j], j) [for x advection]

      or

          (i, z.postpointmesh[k,i,j]) [for v advection]

      but instead just gives that value z.postpointmesh[k,i,j] of the index
      that changes. Hence, the remap procedure must be encoded with an
      understanding of which index is updated in order to complete the push.

      It then makes sense to speak in terms such as:

          z.postpointmesh[0,i,j] + 1 is a contiguous gridpoint to
                       the index z.postpointmesh[0,i,j] in the positive
                       direction

          z.postpointmesh[0,i,j] - 1 is a contiguous gridpoint to
                       the index z.postpointmesh[0,i,j] in the negative
                       direction

      since we are speaking of one index, not a duple as concerns the postpointmesh
      object
    """
    z.postpointmesh[0,:,:] = z.prepointmesh + z.CFL.int
    z.postpointmesh[1,:,:] = np.sign(z.CFL.frac).astype(int) + z.postpointmesh[0,:,:]

    return z.postpointmesh

def remap_step(
        sim_params,
        f_old,
        n,
        z,
        vz
        ):
    """Orchestrates remapping of all advected moving cells to the grid,
    the flowchart of methods being used is the following:

    f_final = convect.remap_step(
               sim_params,
                f_old,
                n,
                z,
                vz
                )

                | |
               \   /
                \ /

    Uf = convect.flux(    ->  c = HOC.correctors(sim_params, z, vz)
        sim_params,
        f_old,           ->   d = derivatives.method(f, z, vz, sim_params)
        z, vz
        )                <-   Uf = sum over q (c[q,:]*d[q,:,:])


                | |
               \   /
                \ /


   remap the appropriate proportion to the nearest neighbor gridpoints
   f_k1 = convect.remap_assignment(
                            f_old,
                            Uf,
                            z.postpointmesh[0,:,:],
                            z,
                            vz,
                            index = 'nearest'
                            )


   remap the remaining proportion to the appropriate contiguous gridpoint
   f_k2 = convect.remap_assignment(
                            f_old,
                            Uf,
                            z.postpointmesh[1,:,:],
                            z,
                            vz,
                            index = 'contiguous'
                            )

   return f_remapped = f_k1 + f_k2

                | |
               \   /
                \ /


   f_final = convect.finalize_density(
                            sim_params,
                            f_remapped,
                            f_final # initialized container of zeros
                            )

    inputs:
    sim_params -- (dict) simulation parameters
    f_old -- (ndarray, dim=2) density from previous time step, full grid
    n  -- (int) current time step
    z -- (instance) phase space variable
    vz -- (instance) generalized velocity variable for phase space variable z

    outputs:
    f_remapped -- (ndarray, dim=2) density with all MCs remapped at final postpoints
                     according to remap rule

    NOTE: we carry assemble the new density as a sum of two
    2D matrices f_k1 and f_k2 which themselves are the associated
    mappings to the postpoints k1 (nearest neighbor grid index to
    exact non-integral postpoint), and k2 (contiguous grid index
    to k1 based on the sign of advection), we use k here for brevity
    but the postpoint pairs are stored in z.postpointmesh[:,i,j] for
    each [i,j]

    The reason a partition is necessary is that a mapping such as

        f_new[k, vz.prepointmesh] += f_old

    does not increment each element, but instead overwrites leaving
    only the final assignment to a particular gridpoint rather than
    adds on top of it. Thus, the separation of k1 and k2 into
    different containers f_k1 and f_k2 permit non-overlapping
    assignments in each individually since every z is pushed by the same
    vz for a given row (column), so overlapping assignments are not
    possible for each index separately. The sum gives the appropriate ]
    total density f_new.
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
        z.postpointmesh[0,:,:],
        z,
        vz,
        index = 'nearest'    # remaps to nearest neighbor index
        )

    # remap to contiguous cell center
    f_k2 = remap_assignment(
        f_old,
        Uf,
        z.postpointmesh[1,:,:],
        z,
        vz,
        index = 'contiguous'    # remaps to contiguous neighbor of above
        )

    f_remapped = f_k1 + f_k2

    # check for density conservationfor every gridpoint

    k = z.postpointmesh

    for i in range(f_old.shape[0]):
        for j in range(f_old.shape[1]):
            if (f_old[i,j] != f_k1[k[0,i,j],j] + f_k2[k[1,i,j], j]):
                print "density not conserved from prepoint [%d, %d], the difference is %g" % (i,j,
                                                                                              f_old[i,j] - (f_k1[k[0,i,j],j] + f_k2[k[1,i,j], j]))


    # global check on density conservation
    DECSKS.lib.density.global_conservation_check(sim_params, f_remapped, n)

    return f_remapped
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

    # evaluate derivatives q = 0, 1, 2, ... N-1 (zeroeth is density itself)

    # calls lib.derivatives.fd or lib.derivatives.fourier based on the
    # HOC specified in etc/params.dat, sim_params['derivative_method']
    # contains the related function handle as a string
    d = eval(sim_params['derivative_method'])(f_old, z, vz, sim_params)

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
    ith nominal flux matrix Uf, shape = (z1.N, z2.N)

    inputs:
    f_old -- (ndarray, dim=2) density from previous time step
    CFL -- (instance) Courant numbers dictating phase space advection
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
        zpostpointmesh,
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

    # TODO this probably has to be switched around for v advection
    # should be index slicing [vz.prepointmesh, zpostpointmesh] for
    # vz = ax, z = vx maaayyyyybbeeee....everything is transposed after all
    if index == 'nearest':
        # mask out negative values
        f_old_ma.mask = mask_neg
        Uf_ma.mask = mask_neg

        f_pos[ zpostpointmesh, vz.prepointmesh ] = f_old_ma - Uf_ma

        # mask out all positive values
        f_old_ma.mask = mask_pos
        Uf_ma.mask = mask_pos

        f_neg[ zpostpointmesh, vz.prepointmesh ] = f_old_ma + Uf_ma

    elif index == 'contiguous':
        # mask out negative values
        f_old_ma.mask = mask_neg
        Uf_ma.mask = mask_neg

        f_pos[ zpostpointmesh, vz.prepointmesh ] = Uf_ma

        # mask out all positive values
        f_old_ma.mask = mask_pos
        Uf_ma.mask = mask_pos

        f_neg[ zpostpointmesh, vz.prepointmesh ] = -Uf_ma

    # "wherever there is negative data, assign f_neg, else assign f_pos
    f_new = np.where(mask_neg == True, f_neg.data, f_pos.data)

    return f_new

def finalize_density(sim_params, f_remapped, f_final, z, vz):
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
    if z.str[0] == 'v':
        DECSKS.lib.domain.velocity_advection_prep(f_final, f_remapped, z, vz) # undo all transpositions
        f_final = np.transpose(f_final)

    # assign all active grid points to grid values on f_final
    f_final[:f_remapped.shape[0], :f_remapped.shape[1]] = f_remapped

    # complete by periodicity at gridpoints [x.N + 1, :],
    # and [:, vx.N + 1] if PBCs on both, then
    f_final[f_remapped.shape[0],:-1] = f_remapped[0,:]
    f_final[:-1, f_remapped.shape[1]] = f_remapped[:,0]

    # complete the diagonal entry [x.N+1, vx.N + 1] with either [0,vx.N+1]
    # or [x.N+1, 0], both will be the same because of dual periodicity.
    f_final[f_remapped.shape[0], f_remapped.shape[1]] = f_final[0, f_remapped.shape[1]]

    if z.str[0] == 'v':
        f_final = np.transpose(f_final)

    return f_final
