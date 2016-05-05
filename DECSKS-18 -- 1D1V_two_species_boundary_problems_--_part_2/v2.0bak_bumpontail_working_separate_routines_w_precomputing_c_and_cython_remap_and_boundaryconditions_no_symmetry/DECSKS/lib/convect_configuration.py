import DECSKS
import numpy as np
import numpy.ma as ma

def scheme(
    f_initial,
    n,
    stage,
    sim_params,
    c,
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
    z -- (instance) phase space variable
    vz -- (instance) generalized velocity for z


    outputs:
    f_final -- (ndarray, ndim = 1 or 2) f[n,:] or f[n,:,:] updated
               after all steps have been completed
    """



    # (0) INITIALIZE FINAL DENSITY CONTAINER AND EXTRACT EVOLVED GRID
    f_final = np.zeros(f_initial.shape)
    f_initial = DECSKS.lib.domain.extract_active_grid(f_initial, z, sim_params)

    # (1) ADVECT DENSITY AND COMPUTE CORRESPONDING FLUXES
    z.postpointmesh = advection_step(stage, z)

    # compute high order fluxes
    Uf = flux(
        sim_params,
        stage,
        c,
        f_initial,
        z, vz
        )

    # (2) APPLY BOUNDARY CONDITIONS AND REMAP DENSITY TO GRID
    f_remapped = remap_step(
                       sim_params,
                       f_initial,
                       Uf,
                       n,
                       stage,
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
def advection_step(stage, z):
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
    z.postpointmesh[0,:,:] = z.prepointmesh + z.CFL.int[stage,:,:]
    z.postpointmesh[1,:,:] = np.sign(z.CFL.numbers[stage,:,:]).astype(int) + z.postpointmesh[0,:,:]

    return z.postpointmesh

def remap_step(
        sim_params,
        f_old,
        Uf_old,
        n,
        stage,
        z,
        vz,
        charge
        ):

    f_copy =  f_old.copy()
    Uf_copy = Uf_old.copy()

    f_copy, Uf_copy, z = \
      eval(sim_params['distribution_function_boundarycondition_orchestrator_handle'][z.str])(
          f_copy, Uf_copy, z, vz, sim_params, charge, k = 0)

    # for now, we do not have a symmetry boundary condition, hence vz.postpointmesh[k,:,:] = vz.prepointmesh
    f_k1 = DECSKS.lib.remap.nearest_gridpoint_assignment(f_copy, Uf_copy, z.postpointmesh[0,:,:], vz.prepointmesh, z.N, vz.N)

    f_old, Uf_old, z = \
      eval(sim_params['distribution_function_boundarycondition_orchestrator_handle'][z.str])(
          f_old, Uf_old, z, vz, sim_params, charge, k = 1)

    # for now, we do not have a symmetry boundary condition, hence vz.postpointmesh[k,:,:] = vz.prepointmesh

    f_k2 = DECSKS.lib.remap.contiguous_gridpoint_assignment(f_old, Uf_old, z.postpointmesh[1,:,:], vz.prepointmesh, z.N, vz.N)    

    f_remapped = f_k1 + f_k2

    return f_remapped

# ........................................................................... #
def flux(
        sim_params,
        stage,
        c,
        f_old,
        z, vz
        ):
    """Computes fluxes Uf for all z.prepointmesh

    inputs:
    sim_params -- (dict) simulation parameters
    c -- (ndarray, ndim=3) array of corrections according to stage s, c[s,:,:] gives
         the 2D set of correctors for stage s, where the indexing is [s,q,j], where
         q labels the qth order correction (q = 0, 1, ... N-1), and j labels the
         normalized velocity, i.e. the CFL = v / v_grid, where v_grid = dx / dt,
         all points i in f[i,j] are characterized with the same CFL label j for
         a given stage s. That is, c[s,:,j] gives the correctors that apply to all i
         in f[i,j]
    s -- (int) current stage of time splitting scheme
    f_old -- (ndarray, dim=2) density from previous time (sub)step
    CFL -- (instance) CFL number
    z -- (instance) phase space variable
    vz -- (instance) generalized velocity corresponding to above phase space variable

    outputs:
    Uf -- (ndarray, dim=2) Normalized fluxes originating
          from every for every [i,j] in z.prepointmesh

           where for each [i,j]:

               Uf[i,j] = sum c[q,j]*d[q,i,j] over q = 0, 1, ... N-1
                       = f_old + (high order corrections)
    """
    # c for configuration variables has been computed in lib.split.scheme
    # and passed as a parameter in this function

    # evaluate derivatives q = 0, 1, 2, ... N-1 (zeroeth is density itself)

    # calls lib.derivatives.fd or lib.derivatives.fourier based on the
    # HOC specified in etc/params.dat, sim_params['derivative_method']
    # contains the related function handle as a string
    d = eval(sim_params['derivative_method'][z.str])(f_old, z, vz, sim_params)

    # compute high order fluxes column-by-column
    Uf = np.zeros(f_old.shape)
    for j in range(vz.N):
        Uf[:,j] = c[:,j].dot(d[:,:,j])

    # enforce flux limiter to ensure positivity and restrict numerical overflow
    Uf = flux_limiter(stage, f_old, Uf, z)

    return Uf

def flux_limiter(
        stage,
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
    Uf_ma = ma.array(z.CFL.frac[stage,:,:])

    # mask negative values, keep positives
    Uf_ma[z.CFL.numbers[stage,:,:] < 0] = ma.masked

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
