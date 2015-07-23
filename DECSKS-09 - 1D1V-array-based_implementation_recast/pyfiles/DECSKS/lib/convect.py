import DECSKS
import numpy as np

def scheme(
    f_old,
    z,n,
    sim_params
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
        z,n,
        sim_params
        )

    return f_new
#---------------------------------------------------------------------------  #
def advection_step(
    f_old,
    z,n,
    sim_params
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

    f_new = remap_step(f_old, CFL, z, n, sim_params)

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
        self.int = np.where(self.numbers >=0, np.floor(self.numbers),
                            np.ceil(self.numbers))

        self.frac = self.numbers - self.int
#---------------------------------------------------------------------------  #
def remap_step(
        f_old,
        CFL,
        z, n,
        sim_params
        ):
    """Orchestrates remapping of all MCs to Eul. mesh by calling RemapRule(*args)

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

    beta = DECSKS.lib.HOC.beta_m(CFL,
                                 sim_params['Bernoulli_numbers'],
                                 sim_params['N'])

    c = np.zeros(sim_params['N']) # if FD, change 'N' to order of FD scheme
    for q in range(1, sim_params['N']):
        c[q] = (-1) ** q * beta[q]

    Uf = flux(
        f_old,
        CFL,
        z, sim_params,
        c
        )

    for i in z.prepoints:
        f_remapped_MC_container = remap_rule(f_old, Uf, CFL, z, sim_params, i, n)
        f_new += f_remapped_MC_container

    if sim_params['BC'].lower() == 'periodic':
        f_new[z.Ngridpoints - 1] = f_new[0]

    return f_new
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
def flux(
        f_old,
        CFL,
        z, sim_params,
        c = None
        ):
    """Computes fluxes Uf for all z.prepoints

    inputs:
    f_old -- (ndarray, dim=1) density from previous time step
    CFL -- (instance) CFL number
    z -- (instance) phase space variable
    sim_params -- (dict) simulation parameters

    outputs:
    Uf -- (ndarray, dim=1) Normalized fluxes for every z.prepoints[i]
    """
    # active gridpoints only
    f_old = f_old[:eval(sim_params['phasespace_vars'][0]).N,
                  :eval(sim_params['phasespace_vars'][1]).N]

    Uf = np.zeros(f_old.shape)
    d = np.zeros([z.N,sim_params['N']])

    # Compute uncorrected fluxes
    Uf = f_old * CFL.frac

    # Compute any high order corrections c1*d1 + c2*d2 + ...
    if sim_params['HOC'] == 'FD':
        d = DECSKS.lib.derivatives.finite_differences_matrix_form(f_old, z, sim_params)
        for q in range(1,sim_params['N']):
            Uf += c[q]*d[:,q] # Flux = G + H.O.C.

    elif sim_params['HOC'] == 'FOURIER':
        if sim_params['WindowedFilter'] == 'YES':
            K = DECSKS.lib.HOC.kernel(z)
        else:
            K = None
        for q in range(1, sim_params['N']):
            d[:,q] = DECSKS.lib.derivatives.trigonometric(f_old,z,q,sim_params, K)
            Uf += c[q]*d[:,q] # Flux = G + H.O.C.
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
