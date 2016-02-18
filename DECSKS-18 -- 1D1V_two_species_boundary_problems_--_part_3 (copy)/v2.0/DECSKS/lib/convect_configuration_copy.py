import DECSKS
import numpy as np
import numpy.ma as ma

def scheme(
    f_initial,
    s,n,
    sim_params, c,
    z,
    vz,
    charge = -1
    ):
    """Solves a collection of 1D advection (in z) equations by convected scheme
    and stacks the results in a 2D matrix

    inputs:
    f_initial -- (ndarray, ndim = 2) f[n-1,:,:] if first
             substep in a splitting algorithm or if none
             else, f[n,:,:]
    n -- (int) time step
    sim_params -- (dict) simulation parameters
    z -- (instance) generic phase space variable (x,y, or physical z)
    vz -- (instance) generalized velocity for z (vx, vy, or physical vz)


    outputs:
    f_final -- (ndarray, ndim = 1 or 2) f[n,:] or f[n,:,:] updated
               after all steps have been completed
    """

    # (0) INITIALIZE FINAL DENSITY CONTAINER AND EXTRACT EVOLVED GRID
    f_final = np.zeros(f_initial.shape)
    f_initial = DECSKS.lib.domain.extract_active_grid(f_initial, sim_params)

    # (1) ADVECT DENSITY AND COMPUTE CORRESPONDING FLUXES
    z.postpointmesh = advection_step(s,z)

    # compute high order fluxes
    Uf = flux(
        sim_params, c,
        s,
        f_initial,
        z, vz
        )

    # (2) APPLY BOUNDARY CONDITIONS AND REMAP DENSITY TO GRID
    f_remapped = remap_step(
                       sim_params,
                       f_initial,
                       Uf,
                       s, n,
                       z,
                       vz,
                       charge
                       )

    # (3) COLLISION STEP (NOT YET IMPLEMENTED)
    # f_new = DECSKS.lib.collisions.collisiontype(f_old, z, n)

    # (4) RETURN FINAL DESTINY (density*)
    f_final = finalize_density_absorbing(sim_params, f_remapped, f_final, z, vz)

    return f_final
#---------------------------------------------------------------------------  #
def advection_step(s,z):
    """Pushes each z.prepointmesh (index) value by the advection *along each
    column j* (i.e. constant vz.prepointvaluemesh)
    as prescribed by its generalized velocity vz.prepointmeshvalues[:,j].

    This is computed in one calculation by:

        vz.prepointmeshvalues * dt / z.width = CFL.numbers

    the integer part, CFL.int, is pushed here. The residual fraction
    CFL.frac is remapped according to lib.convect.remap_step in step (2) in
    the orchestrator routine, lib.convect.scheme, above.

    inputs:
    s -- (int) stage of the split scheme
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


    USAGE NOTE: the map z.postpointmesh[k,i,j] for k = 0 or 1 gives the postpoint
      of the advecting index, it does not give the relevant postpoint duple
      in phase space:

          (z.postpointmesh[k,i,j], j) [for x advection]

      or

          (i, z.postpointmesh[k,i,j]) [for v advection]

      but instead just gives that value z.postpointmesh[k,i,j] of the index
      that changes. Hence, the remap procedure must be encoded with an
      understanding of which index is updated in order to complete the push.

      It then makes sense (broadcasting is clear) to speak in terms such as:

          z.postpointmesh[0,i,j] + 1 is a contiguous gridpoint to
                       the index z.postpointmesh[0,i,j] in the positive
                       direction

          z.postpointmesh[0,i,j] - 1 is a contiguous gridpoint to
                       the index z.postpointmesh[0,i,j] in the negative
                       direction

      since we are speaking of one index, not a duple as concerns the postpointmesh
      object
    """
    # NOTE: we must use the sign on the CFL numbers as it is not interchangeable with the
    # velocity values given high order splitting methods require negative splitting coefficients
    # hence the direction of the trajectory is not just dependent on vz.prepointvaluemesh but
    # also on vz.prepointvaluemesh*(split_coeff*t.width)


    z.postpointmesh[0,:,:] = z.prepointmesh + z.CFL.int[s,:,:]
    z.postpointmesh[1,:,:] = np.sign(z.CFL.numbers[s,:,:]).astype(int) + z.postpointmesh[0,:,:]

    return z.postpointmesh

def remap_step(
        sim_params,
        f_template,
        Uf_template,
        s, n,
        z,
        vz,
        charge
        ):
    """Orchestrates remapping of all advected moving cells to the grid,

    First, we map the appropriate proportion to the nearest grid point (k1) to
    the exact (non-integral in general) postpoint.

    Then, we map the remaining fraction to the contiguous grid point (k2)
    defined to be in the same direction of travel as the advection

    inputs:
    sim_params -- (dict) simulation parameters
    f_template = f_old --
              (ndarray, dim=2) density from previous time step, full grid.
              We call this a template because it is used as a template
              in using it for achieving the boundary conditions.

              (1) we copy a template and send it to the boundary conditions
                  when preparing to remap density to the k1 postpoints.
                  To achieve the BCs, we modify the density and fluxes
                  (e.g. zeroing out any absorbed at walls), then remap
                  to a container f_k1

              (2) we repeat the same as above, but this time in preparation
                  for mapping at the contiguous grid point k2. If we used
                  the same density as was returned in (1), it already has
                  been mangled to acheive boundary conditions. Thus, we
                  need a fresh template. Since there are only two postpoints
                  (k1 and k2), we do not need a copy of the original template
                  but can instead just modify the template as its utility
                  ends here.

    Uf_template = Uf_old -- (ndarrray, ndim=2), the flux values for the
                            full grid. See the description for f_template
                            on why this is named as Uf_template rather
                            than Uf_old.

    s -- (int) current splitting stage; used to access the pre-computed
            correctors c, which depend on the stage. Note, for physical
            velocity (e.g. in lib.convect_velocity) we must compute the
            correctors in each pass as these are not tied to any grid
            hence are difference from one full time step to the next
    n  -- (int) current time step
    z -- (instance) phase space variable
    vz -- (instance) generalized velocity variable for phase space variable z
    charge -- (int) -1 or +1, indicates charge species

    outputs:
    f_remapped -- (ndarray, dim=2) density with all densities reallocated
                 according to their postpoint maps:

      postpoint map = (z.postpointmesh[k1 or k2, :,:], vz.prepointmesh)

    here, we emphasize that it is z that is being advected, hence
    these have distinct postpoints. While we may modify the values
    of the generalized velocity in special circumstances (e.g.
    symmetry boundary condition), the velocities of the prepoints
    are interpreted physically as being the same before (prepoints)
    and after (postpoints) since we are evolving the system with
    splitting routine (one variable evolved at a time)
    """

    f_copy = np.copy(f_template)
    Uf_copy = np.copy(Uf_template)

    # Prepare for remapping to k1

    # apply boundary conditions for density packets reaching [k1, j]
    # here, we pass k = 0 as a parameter, as it refers to z.postpointmesh[k,:,:]
    # and k1 corresponds to the storage k = 0
    f_copy, Uf_copy = \
      eval(sim_params['boundarycondition_function_handle'][z.str])(
          f_copy, Uf_copy, z, vz, sim_params, charge, k = 0)

    Uf_nonneg, Uf_neg = DECSKS.lib.remap.sift(Uf_copy, z.CFL.numbers[s,:,:])

    # remap all [i,j] to postpoints [k1, j], we assign per the piecewise rule:
    #
    #        f_k1[k1,j] =  f_copy[i,j] - Uf_copy[i,j] if CFL >= 0
    #
    #                      f_copy[i,j] + Uf_copy[i,j] if CFL < 0
    #
    # we accomplish the above through the following set of operations in order to minimize the computational cost
    f_k1 = np.zeros_like(f_template)
    f_k1  = DECSKS.lib.remap.assignment(f_copy, z.postpointmesh[0,:,:], vz.postpointmesh[0,:,:], f_k1.shape[0], f_k1.shape[1])
    f_k1 += DECSKS.lib.remap.assignment(Uf_neg, z.postpointmesh[0,:,:], vz.postpointmesh[0,:,:], f_k1.shape[0], f_k1.shape[1])
    f_k1 -= DECSKS.lib.remap.assignment(Uf_nonneg, z.postpointmesh[0,:,:], vz.postpointmesh[0,:,:], f_k1.shape[0], f_k1.shape[1])

    print f_k1
    # store an array of indices which indicate the "special" entries, i.e. those that are around the edges and whcih require both
    # the partner flux and the non-partner fluxes

    # edge_postpoints = np.where(0 <= z.postpointmesh[0,:,:] <= 1)
    # gives the (i,j) pairs, only need to base diagnostic off of one k1 or k2


    # Prepare for remapping to k2
    # We do not need the information in f_template, and Uf_template hereafter, so we may modify these directly

    # apply boundary conditions for density packets reaching [k2, j]
    # here, we pass k = 0 as a parameter, as it refers to z.postpointmesh[k,:,:]
    # and k2 corresponds to the storage k = 1
    f_template, Uf_template = \
      eval(sim_params['boundarycondition_function_handle'][z.str])(
          f_template, Uf_template, z, vz, sim_params, charge, k = 1)

    Uf_nonneg, Uf_neg = DECSKS.lib.remap.sift(Uf_template, z.CFL.numbers[s,:,:])

    # remap all [i,j] to postpoints [k2, j], we assign per the piecewise rule:
    #
    #        f_k2[k2,j] =  -Uf_template[i,j] if CFL >= 0
    #
    #                      +Uf_template[i,j] if CFL < 0
    #
    # we accomplish the above through the following set of operations in order to minimize the computational cost
    f_k2 = np.zeros_like(f_template)
    f_k2 -= DECSKS.lib.remap.assignment(Uf_neg, z.postpointmesh[1,:,:], vz.postpointmesh[1,:,:], f_k2.shape[0], f_k2.shape[1])
    f_k2 += DECSKS.lib.remap.assignment(Uf_nonneg, z.postpointmesh[1,:,:], vz.postpointmesh[1,:,:], f_k2.shape[0], f_k2.shape[1])

    f_remapped = f_k1 + f_k2

    return f_remapped

# ........................................................................... #
def flux(
        sim_params, c,
        s,
        f_old,
        z, vz
        ):
    """Computes fluxes Uf for all z.prepointmesh

    inputs:
    sim_params -- (dict) simulation parameters
    f_old -- (ndarray, dim=2) density from previous time (sub)step
    z -- (instance) phase space variable
    vz -- (instance) generalized velocity variable

    outputs:
    Uf -- (ndarray, dim=2) Normalized fluxes originating
          from every for every [i,j] in z.prepointmesh

           where for each [i,j]:

               Uf[i,j] = sum c[q,j]*d[q,i,j] over q = 0, 1, ... N-1
                       = f_old + (high order corrections)
    """

    # the correctors c[s,q,j] were calculated prior to the simulation start
    # inside lib.split

    # evaluate derivatives q = 0, 1, 2, ... N-1 (zeroeth is density itself)

    # calls lib.derivatives.fd or lib.derivatives.fourier based on the
    # HOC specified in etc/params.dat, sim_params['derivative_method']
    # contains the related function handle as a string
    d = eval(sim_params['derivative_method'][z.str])(f_old, z, vz, sim_params)

    # compute high order fluxes column-by-column
    Uf = np.zeros(f_old.shape)
    for j in range(vz.N):
        Uf[:,j] = c[s,:,j].dot(d[:,:,j])

    # enforce flux limiter to ensure positivity and restrict numerical overflow
    Uf = flux_limiter(s, f_old, Uf, z)

    return Uf

def flux_limiter(
        s,
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
    Uf_ma = ma.array(z.CFL.frac[s,:,:])

    # mask negative values, keep positives
    Uf_ma[z.CFL.frac[s,:,:] < 0] = ma.masked

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

def finalize_density_periodic(sim_params, f_remapped, f_final, z, vz):
    """
    Returns a final density assuming x and vx have periodic boundaries.
    Recall that (z1.N, z2.N) moving cells are evolved.

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

    # assign all active grid points to grid values on f_final, the indexing is unnecessary
    # but included for transparency, f_final (total grid) contains an extra column and
    # row as compared to f_remapped (active grid)
    f_final[:f_remapped.shape[0], :f_remapped.shape[1]] = f_remapped

    # complete by periodicity at gridpoints [x.N + 1, :],
    # and [:, vx.N + 1] if PBCs on both, then
    f_final[f_remapped.shape[0],:-1] = f_remapped[0,:]
    f_final[:-1, f_remapped.shape[1]] = f_remapped[:,0]

    # complete the diagonal entry [x.N+1, vx.N + 1] with either [0,vx.N+1]
    # or [x.N+1, 0], both will be the same because of dual periodicity.
    f_final[f_remapped.shape[0], f_remapped.shape[1]] = f_final[0, f_remapped.shape[1]]

    return f_final

def finalize_density_absorbing(sim_params, f_remapped, f_final, z, vz):
    """
    finalizes density function assuming x nonperiodic, vx periodic
    """

    # assign all active grid points to grid values on f_final
    f_final[:f_remapped.shape[0], :f_remapped.shape[1]] = f_remapped

    # complete vx by periodicity at gridpoints [:, vx.N + 1], indexing f_remapped.shape[1] is the final
    # gridpoint given zero based indexing.
    f_final[:, f_remapped.shape[1]] = f_remapped[:,0]

    return f_final
