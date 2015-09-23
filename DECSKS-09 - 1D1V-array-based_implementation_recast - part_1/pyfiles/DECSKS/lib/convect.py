import DECSKS
import numpy as np
import numpy.ma as ma

def scheme(
    f_old,
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
    f_new -- (ndarray, ndim = 1 or 2) f[n,:] or f[n,:,:] updated
    """

    f_new = np.zeros(f_old.shape)
    f_new = advection_step(
        f_old,
        n,
        sim_params,
        z,
        vz
        )

    return f_new
#---------------------------------------------------------------------------  #
def advection_step(
    f_old,
    n,
    sim_params,
    z,
    vz
    ):
    """Pushes phase space z.MCs by an integral number of cells (CFL.int)
    then remaps MCs according to RemapStep(*args)

    inputs:
    f_old -- (ndarray, dim=1) f(z1,z2=const,z3=const,...), i.e. we split
        the PDE as 1D advection equations
    z -- (instance) phase space variable
    n -- (int) current time step
    sim_params -- (dict) simulation parameters

    outputs:
    f_new -- (ndarray, dim=1) density after time step has been taken completely
    """
    # advection step
    CFL = CourantNumber(z)  # instance with attr. CFL.int, CFL.frac
    z.postpointmesh = z.prepointmesh + CFL.int

    # apply boundary conditions
    z.postpointmesh = DECSKS.lib.boundaryconditions.periodic(z)

    # collision step (Not yet implemented), e.g. collisiontype = Coulomb, ...
    # f_new = DECSKS.lib.collisions.collisiontype(f_old, z, n)

    f_new = remap_step(f_old,
                       CFL,
                       n,
                       sim_params,
                       z,
                       vz
                       )

    return f_new
#---------------------------------------------------------------------------  #
class CourantNumber:
    """Returns a CFL number instance of the phasespace variable z

    inputs:
    z -- (instance) phase space variable from class lib.domain.Setup

    outputs:
    self -- (instance) CFL number ascribed to variable z convection

    Note: broadcasting ensures self.numbers is the same shape as z.MCs
    """
    def __init__(self, z):

        self.numbers = (z.MCs - z.prepointvaluemesh) / z.width

        # if >= 0 , self.int = floor(self.numbers), else ceil(self.numbers)
        # i.e. the sign of CFL.frac agrees with the sign of vz
        self.int = np.where(self.numbers >=0, np.floor(self.numbers),
                            np.ceil(self.numbers))

        self.frac = self.numbers - self.int
#---------------------------------------------------------------------------  #
def remap_step(
        f_old,
        CFL,
        n,
        sim_params,
        z,
        vz
        ):
    """Orchestrates remapping of all MCs to Eul. mesh by calling remap_rule(*args)

    inputs:
    f_old -- (ndarray, dim=1) density from previous time step
    Uf -- (ndarray, dim=1) the normalized flux used in CS update
    z -- (instance) phase space variable
    n  -- (int) current time step

    outputs:
    f_new -- (ndarray, dim=1) density with all MCs remapped at final z.postpoints
        final postpoints

    note: z.N MCs are evolved, for peroidic BCs, z.N = len(f_old) - 1
                               and we enforce the BC: f[z.Ngridpoints-1] = f[0])
                               at the conclusion of the timestep (this function)
          for all other boundary problems, z.N = z.Ngridpoints = len(f_old)

          hence, f_old[:z.N] is passed instead of f_old

    """
    f_new = np.zeros(f_old.shape)

    Uf = flux(
        sim_params,
        f_old,
        CFL,
        z, vz
        )

    for j in vz.prepoints:
        f_remapped_MC_containers = remap_rule_2D(f_old,
                                              n,
                                              CFL,
                                              Uf,
                                              sim_params,
                                              z,
                                              vz
                                              )

    # TODO remove the following in favor of the above
    # TODO need to have a matrix of k values,


    for i in z.prepoints:
        f_remapped_MC_container = remap_rule(f_old, Uf, CFL, z, sim_params, i, n)
        f_new += f_remapped_MC_container


    # TODO needs to select based on the active grid direction
    if sim_params['BC'].lower() == 'periodic':
        f_new[z.Ngridpoints - 1] = f_new[0]

    return f_new
# ........................................................................... #
def flux(
        sim_params,
        f_old,
        CFL,
        z,
        vz
        ):
    """Computes fluxes Uf for all z.prepoints

    inputs:
    f_old -- (ndarray, dim=2) density from previous time (sub)step
    CFL -- (instance) CFL number
    z -- (instance) phase space variable
    sim_params -- (dict) simulation parameters

    outputs:
    Uf -- (ndarray, dim=1) Normalized fluxes for every z.prepoints[i]
    """
    # active gridpoints only, note for all but PBCs this is the entire grid
    # so no slicing is actually happening
    f_old = f_old[:eval(sim_params['phasespace_vars'][0]).N,
                  :eval(sim_params['phasespace_vars'][1]).N]

    # TODO need to make this general, so that separate assignments for Bx and Bv
    c = DECSKS.lib.HOC.correctors(sim_params, CFL, z)
    d = np.zeros([sim_params['N'], eval(sim_params['phasespace_vars'][0]).N,
                  eval(sim_params['phasespace_vars'][1]).N])


    Uf = np.zeros(f_old.shape)

    # assemble string name of chosen derivative function

    # call routine (fd or fourier) to compute all derivatives
    d[0,:,:] = f_old
    d = eval(sim_params['derivative_method'])(d, f_old, z, sim_params)

    for j in range(vz.N):
        Uf[:,j] = c[:,j].dot(d[:,:,j])

    # flux limiter to ensure positivity
    Uf = flux_limiter_2D(
        f_old,
        CFL,
        Uf
        )

    return Uf
# ........................................................................... #
def flux_limiter(
        i,
        f_old,
        CFL,
        Uf):
    """Applies the flux (numerical and positivity) limiter to the
    ith nominal flux Uf[i]

    inputs:
    i -- (int) MC prepoint
    f_old -- (ndarray, dim=1) density from previous time step
    CFL -- (instance) Courant numbers characterizing MC convection
    Uf -- (ndarray, dim=1) the normalized flux (high order or not) used in the CS update

    outputs:
    Uf[i] -- (float) final normalized flux for MC originating at z.prepoints[i]
            after numerical and positivity limiter has been applied

    """
    if CFL.frac[i] >= 0:
        Uf[i] = min( [ max([0,Uf[i]]), 1.0*f_old[i]] )
    else: # U < 0
        Uf[i] = min( [ max([-1.0*f_old[i], Uf]), 0 ] )

    return Uf[i]
# ........................................................................... #
def flux_limiter_2D(
        f_old,
        CFL,
        Uf):
    """Applies the flux (numerical and positivity) limiter to the
    ith nominal flux matrix Uf = (z1.N, z2.N)

    inputs:
    f_old -- (ndarray, dim=2) density from previous time step
    CFL -- (instance) Courant numbers characterizing MC convection
    Uf -- (ndarray, dim=2) the normalized flux (high order or not) used in the CS update

    outputs:
    Uf -- (ndarray, dim=2) final normalized flux for MC originating at prepoints[i,j]
            after numerical and positivity limiter has been applied

    Note: the masked implementation here takes about 285 ms for a 768 x 1536 grid
          the looped implementation over all i,j takes 3.85 s for the whole grid

    """
    # local masked array (ma) copy of CFL.frac used to leave CFL.frac unchanged
    Uf_ma = ma.array(CFL.frac)

    # mask negative values, keep positives
    Uf_ma[CFL.frac < 0] = ma.masked

    # assign the mask to a zero matrix of same size for pairwise ma.minimum/maximum operations below
    zeros = ma.zeros(Uf_ma.shape)
    zeros.mask = Uf_ma.mask

    # operating only on positive values
    Uf_pos = ma.minimum(ma.maximum(zeros, Uf), 1.0*f_old)

    # mask positives, keep negatives
    Uf_ma.mask = np.logical_not(Uf_ma.mask)
    zeros.mask = Uf_ma.mask

    #operating only on negative values
    Uf_neg = ma.minimum(ma.maximum(-1.0*f_old, Uf), zeros)

    # combine masked values in a single matrix Uf_final
    Uf = np.where(Uf_neg.mask == False, Uf_neg.data, Uf_pos.data)

    return Uf
# ........................................................................... #
def remap_rule_2D(
        f_old,
        n,
        CFL,
        Uf,
        sim_params,
        z,
        vz
        ):
    """Remaps the iz_th MC to Eulerian mesh at indices k1 and k2

    inputs:
    f_old -- (ndarray, dim=1) density from previous time step
    CFL -- (instance) Courant numbers characterizing MC convection
    Uf -- (ndarray, dim=1) the normalized flux used in the CS update
    z -- (instance) phase space variable whose MCs are convected
    i  -- (int) MC prepoint index
    n -- (int) current time step

    outputs:
    f_remapped_MC_container -- (ndarray, dim=1) density container with remapped MCs at final
                                 postpoints
    """
    # active gridpoints only
    f_old = f_old[0:eval(sim_params['phasespace_vars'][0]).N,
                  0:eval(sim_params['phasespace_vars'][1]).N]

    Uf[i] = flux_limiter(i, f_old, CFL, Uf)

    f_remapped_MC_container = np.zeros(z.Ngridpoints)
    k1, k2 = DECSKS.lib.boundaryconditions.periodic_old(z, i, Uf)
    # remap assignment
    if Uf[i] > 0:
        f_remapped_MC_container[k1] = f_old[ z.prepoints[i] ] - Uf[ z.prepoints[i] ]
        f_remapped_MC_container[k2] = Uf[z.prepoints[i]]

    elif Uf[i] < 0:
        f_remapped_MC_container[k1] = f_old[z.prepoints[i]] + Uf[z.prepoints[i]]
        f_remapped_MC_container[k2] = -Uf[z.prepoints[i]]

    else: #
        f_remapped_MC_container[z.postpoints[i]] = f_old[z.prepoints[i]]

    DECSKS.lib.density.conservation_check(f_remapped_MC_container, f_old, i, n)

    return f_remapped_MC_container
# ........................................................................... #
def remap_rule(
        f_old,
        Uf,
        CFL,
        z, sim_params,
        i,n,
        ):
    """Remaps the iz_th MC to Eulerian mesh at indices k1 and k2

    inputs:
    f_old -- (ndarray, dim=1) density from previous time step
    CFL -- (instance) Courant numbers characterizing MC convection
    Uf -- (ndarray, dim=1) the normalized flux used in the CS update
    z -- (instance) phase space variable whose MCs are convected
    i  -- (int) MC prepoint index
    n -- (int) current time step

    outputs:
    f_remapped_MC_container -- (ndarray, dim=1) density container with remapped MCs at final
                                 postpoints
    """
    # active gridpoints only
    f_old = f_old[0:eval(sim_params['phasespace_vars'][0]).N,
                  0:eval(sim_params['phasespace_vars'][1]).N]

    f_remapped_MC_container = np.zeros(z.Ngridpoints)
    k1, k2 = DECSKS.lib.boundaryconditions.periodic_old(z, i, Uf)
    # remap assignment
    if Uf[i] > 0:
        f_remapped_MC_container[k1] = f_old[ z.prepoints[i] ] - Uf[ z.prepoints[i] ]
        f_remapped_MC_container[k2] = Uf[z.prepoints[i]]

    elif Uf[i] < 0:
        f_remapped_MC_container[k1] = f_old[z.prepoints[i]] + Uf[z.prepoints[i]]
        f_remapped_MC_container[k2] = -Uf[z.prepoints[i]]

    else: #
        f_remapped_MC_container[z.postpoints[i]] = f_old[z.prepoints[i]]

    DECSKS.lib.density.conservation_check(f_remapped_MC_container, f_old, i, n)

    return f_remapped_MC_container
