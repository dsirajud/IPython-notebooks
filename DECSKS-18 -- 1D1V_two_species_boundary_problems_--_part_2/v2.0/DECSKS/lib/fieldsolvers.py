import numpy as np
import numpy.linalg as LA
import DECSKS

def compute_electric_field_fd_nonperiodic(fe, fi, x, vx, n, sim_params):

    phi = Poisson_DBC_1D1V_2S(fe, fi, x, vx, n, sim_params)

    # currently the finite difference weight matrix W_dn1 is a 6th order LTE to
    # match the 6th order LTE on the Poisson solve
    dphi = 1 / x.width ** 1 * sim_params['W_dn1_LTE6'].dot(phi)
    Ex = -dphi

    return Ex

def compute_electric_field_fd_dirichlet(fe, fi, x, vx, n, sim_params):

    phi = Poisson_DBC_1D1V_2S(fe, fi, x, vx, n, sim_params)

    # currently the finite difference weight matrix W_dn1 is a 6th order LTE to
    # match the 6th order LTE on the Poisson solve
    dphi = 1 / x.width ** 1 * sim_params['W_dn1_LTE6'].dot(phi)
    Ex = -dphi

    return Ex

def compute_electric_field_fd_periodic(fe, fi, x, vx, n, sim_params):

    phi = Poisson_PBC_1D1V_2S(fe, fi, x, vx, n, sim_params)

    # currently the finite difference weight matrix W_dn1 is a 6th order LTE to
    # match the 6th order LTE on the Poisson solve
    dphi = 1 / x.width ** 1 * sim_params['W_dn1_LTE6'].dot(phi)
    Ex = -dphi

    return Ex

def compute_electric_field_fourier_periodic(fe, fi, x, vx, n, sim_params):
    """Computes self-consistent electric field E by solving Gauss' law
    using FFT/IFFT for two species.

    inputs:
    ni -- (float) uniform background density of ions,
                  in the future can take an input fi, to compute ni
    f -- (ndarray, dim=2) electron density fe(x,v,n) at time step t^n
                          used to compute ne(x,n) at time step t^n
    x -- (instance) spatial variable
    vx -- (instance) velocity variable
    n -- (int) time step number, t^n
    sim_params -- (dict) simulation parameters dictionary
        sim_params['xi'] -- (ndarray, ndim=1) wave number vector over x

    outputs:
    E -- (ndarray,dim=2) electric field, E(x,v) = E(x) at time t^n for all (i,j)
    """
    fe = DECSKS.lib.domain.extract_active_grid(fe[n,:,:], x, sim_params)
    fi = DECSKS.lib.domain.extract_active_grid(fi[n,:,:], x, sim_params)

    n_total = single_integration(fi - fe, of = x, wrt = vx)

    Fn_total = np.fft.fft(n_total) # transformed density
    FE = np.zeros(n_total.shape, dtype = complex)

    FE[1:] = 1 / (1j * sim_params['xi']['x'][1:]) * Fn_total[1:]
    # implicit here is that FE[0] = 0, i.e. we solved for only the fluctuating
    # portion of the field. If the boundary conditions are *not* periodic
    # then some adjustment will need to be factored in (cf. notebook
    # DECSKS-05 for discussion on how to do this)

    E = np.real(np.fft.ifft(FE))

    # extend for all [i,j]
    E = np.outer(E, np.ones([1, vx.N]))

    return E

def Gauss1D1V(ni, f, x, vx, n, sim_params):
    """Computes self-consistent electric field E by solving Gauss' law
    using FFT/IFFT.

    ** This is the same routine as Gauss, but it makes
       vx.N copies of the electric field vector. This routine
       could be used to calculate electric field related quantities
       such as electrostatic energy, but we would need to only compute
       over on column in order to not multiply the quantities artificially
       by a factor of vx.N

    inputs:
    ni -- (float) uniform background density of ions,
                  in the future can take an input fi, to compute ni
    f -- (ndarray, dim=2) electron density fe(x,v,n) at time step t^n
                          used to compute ne(x,n) at time step t^n
    x -- (instance) spatial variable
    vx -- (instance) velocity variable
    n -- (int) time step number, t^n
    sim_params -- (dict) simulation parameters dictionary
        sim_params['xi'] -- (ndarray, ndim=1) wave number vector over x

    outputs:
    E -- (ndarray,dim=2) electric field, E(x,v) = E(x) at time t^n for all (i,j)
    """
    f = DECSKS.lib.domain.extract_active_grid(f[n,:,:], x, sim_params)
    ne = single_integration(f, of = x, wrt = vx)
    n_total = ni - ne

    Fn_total = np.fft.fft(n_total) # transformed density
    FE = np.zeros(ne.shape, dtype = complex)

    FE[1:] = 1 / (1j * sim_params['xi']['x'][1:]) * Fn_total[1:]
    # implicit here is that FE[0] = 0, i.e. we solved for only the fluctuating
    # portion of the field. If the boundary conditions are *not* periodic
    # then some adjustment will need to be factored in (cf. notebook
    # DECSKS-05 for discussion on how to do this)

    E = np.real(np.fft.ifft(FE))

    # extend for all [i,j]
    E = np.outer(E, np.ones([1, vx.N]))

    return E

def Poisson_PBC(ni, f,
                x, vx, n,
                sim_params):
    """6th order LTE finite difference Poisson solver for periodic BCs

    inputs:
    ni -- (float) uniform background density of ions,
                  in the future can take an input fi, to compute ni
    f -- (ndarray, dim=2) electron density fe(x,v,n) at time step t^n
                          used to compute ne(x,n) at time step t^n
    x -- (instance) spatial variable
    v -- (instance) velocity variable
    n -- (int) time step number, t^n


    outputs:
    phi -- (ndarray,dim=1) scalar potential, phi(x) at time t^n,
           for i = 0, 1, ... , x.N - 1, one full period
    """
    f = DECSKS.lib.domain.extract_active_grid(f[n,:,:], x, sim_params)

    # charge densities
    ne = single_integration(f, of = x, wrt = vx)
    n_total = ne - ni

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = B*n + phi_BC

    # label the RHS as b = dx ** 2 * B*n
    b = x.width ** 2 * sim_params['Poisson_6th_order_FD_solver_matrices']['B'].dot(n_total)

    # solve D*phi = b
    phi = LA.solve(sim_params['Poisson_6th_order_FD_solver_matrices']['D'], b)

    # PBCs do not produce unique solutions but a family of solutions with arbitrary integration constant
    # that corresponds to a DC offset phi_avg, recenter so that phi_avg = 0

    phi_avg = np.sum(phi) * x.width / x.L
    phi -= phi_avg

    return phi

def Poisson_PBC_6th_1D1V(ni, f,
                x, vx, n,
                sim_params):
    """6th order LTE finite difference Poisson solver for periodic BCs

    inputs:
    ni -- (float) uniform background density of ions,
                  in the future can take an input fi, to compute ni
    f -- (ndarray, dim=2) electron density fe(x,v,n) at time step t^n
                          used to compute ne(x,n) at time step t^n
    x -- (instance) spatial variable
    v -- (instance) velocity variable
    n -- (int) time step number, t^n


    outputs:
    phi -- (ndarray,dim=2) scalar potential, phi(x,v) = phi(x) at time t^n,
           for i = 0, 1, ... , x.N - 1, one full period
    """
    f = DECSKS.lib.domain.extract_active_grid(f[n,:,:], x, sim_params)

    # charge densities
    ne = single_integration(f, of = x, wrt = vx)
    n_total = ne - ni

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = B*n + phi_BC

    # label the RHS as b = dx ** 2 * B*n
    b = x.width ** 2 * sim_params['Poisson_6th_order_FD_solver_matrices']['B'].dot(n_total)

    # solve D*phi = b
    phi = LA.solve(sim_params['Poisson_6th_order_FD_solver_matrices']['D'], b)

    # PBCs do not produce unique solutions but a family of solutions with arbitrary integration constant
    # that corresponds to a DC offset phi_avg, recenter so that phi_avg = 0
    phi_avg = np.sum(phi) * x.width / x.L
    phi -= phi_avg

    # generate the 2D map for every [i,j], note that each row is constant
    phi = np.outer(phi, np.ones([1,vx.N]))

    return phi

def Poisson_PBC_1D1V_2S(fe, fi,
                x, vx, n,
                sim_params):
    """6th order LTE finite difference Poisson solver for periodic BCs

    The signature 2S = "two species"

    inputs:
    fe -- (ndarray, dim=2) electron density fe(x,v,n) at time step t^n
                          used to compute ne(x,n) at time step t^n
    fi -- (ndarray, dim=2) ion density fe(x,v,n) at time step t^n
                          used to compute ni(x,n) at time step t^n
    x -- (instance) spatial variable
    v -- (instance) velocity variable
    n -- (int) time step number, t^n


    outputs:
    phi -- (ndarray,dim=2) scalar potential, phi(x,v) = phi(x) at time t^n,
           for i = 0, 1, ... , x.N - 1, one full period
    """
    fe = DECSKS.lib.domain.extract_active_grid(fe[n,:,:], x, sim_params)
    fi = DECSKS.lib.domain.extract_active_grid(fi[n,:,:], x, sim_params)

    # Poisson eq. has -(charge density) = ne - ni
    n_total = single_integration(fe - fi, of = x, wrt = vx)

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = B*n + phi_BC

    # label the RHS as b = dx ** 2 * B*n
    b = x.width ** 2 * sim_params['Poisson_6th_order_FD_solver_matrices']['B'].dot(n_total)

    # solve D*phi = b
    phi = LA.solve(sim_params['Poisson_6th_order_FD_solver_matrices']['D'], b)

    # boundary conditions: periodic

    # periodic do not produce unique solutions but a family of solutions with arbitrary integration constant
    # that corresponds to a DC offset phi_avg, recenter so that phi_avg = 0

    phi_avg = np.sum(phi) * x.width / x.L

    phi -= phi_avg
    phi_avg = np.sum(phi) * x.width / x.L

    # generate the 2D map for every [i,j], note that each row is constant
    phi = np.outer(phi, np.ones([1,vx.N]))

    return phi

def Poisson_DBC_1D1V_2S(fe, fi,
                x, vx, n,
                sim_params):
    """6th order LTE finite difference Poisson solver for Dirichlet BCs

    see notebook DECSKS-04 for details on construction

    https://github.com/dsirajud/IPython-notebooks/blob/master/
    DECSKS-04%20--%20Design%20of%20a%206th%20order%20FD%20Poisson%20solver/
    DECSKS-04%20--%20Design%20of%20a%206th%20order%20FD%20Poisson%20solver.ipynb

    The signature 2S = "two species"

    inputs:
    fe -- (ndarray, dim=2) electron density fe(x,v,n) at time step t^n
                          used to compute ne(x,n) at time step t^n
    fi -- (ndarray, dim=2) ion density fe(x,v,n) at time step t^n
                          used to compute ni(x,n) at time step t^n
    x -- (instance) spatial variable
    v -- (instance) velocity variable
    n -- (int) time step number, t^n


    outputs:
    phi -- (ndarray,dim=2) scalar potential, phi(x,v) = phi(x) at time t^n,
           for i = 0, 1, ... , x.N - 1, one full period
    """
    fe = DECSKS.lib.domain.extract_active_grid(fe[n,:,:], x, sim_params)
    fi = DECSKS.lib.domain.extract_active_grid(fi[n,:,:], x, sim_params)

    # Poisson eq. has -(charge density) = ne - ni
    n_total = single_integration(fe - fi, of = x, wrt = vx)

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = b
    #
    # where     b = x.width ** 2 * B*n + phi_BC

    # boundary conditions: dirichlet

    phi_DBC = np.zeros(x.N)

    print sim_params['BC']['x']['phi']['lower']
    print sim_params['BC']['x']['phi']['upper']
    phi_DBC[0] = sim_params['BC']['x']['phi']['lower']
    phi_DBC[-1] = sim_params['BC']['x']['phi']['upper']

    b = x.width ** 2 * sim_params['Poisson_6th_order_FD_solver_matrices']['B'].dot(n_total) + phi_DBC

    # solve D*phi = b
    phi = LA.solve(sim_params['Poisson_6th_order_FD_solver_matrices']['D'], b)

    # generate the 2D map for every [i,j], note that each row is constant
    phi = np.outer(phi, np.ones([1,vx.N]))

    return phi

# INSERT SOLVERS for NBC/DBC, DBC/NBC, and NBC-DBC (Cauchy BC at left boundary <- NBC/NBC)
# need to set up Poisson matrices accordingly in lib.read and access per
# sim_params['Poisson_6th_order_FD_solver_matrices']['B'] and
# sim_params['Poisson_6th_order_FD_solver_matrices']['D']





















def single_integration(f, of = None, wrt = None):
    """integrates once a two variable
    function, i.e. computes integral f(z,wrt) d(wrt) = f(z),
    or integral f(wrt) d(wrt) = F. the keyword
    'of' is the unintegrated variable such that f is a function
    'of' that variable after it is integrated with respect
    to the variable 'wrt'. Momentarily writing of = z, the
    returned integrated function would
    then be F = F(z). If of = None, then the return is
    F = sum(F)*wrt.width = constant.

    Note: there is no need for conditional checks here,
    if we wish to integrate along rows we specificy axis = 0
    in numpy.sum. If the density passed is 1D, axis = 0 still
    adds up all the entries as needed. The axis argument
    indicates the 'direction' of summing.

    for a 2D matrix = (rows, cols):
    axis = 0 would sum every row for each column
    axis = 1 would sum every column for each row

    inputs:
    f -- (ndarray, ndim = 1,2) density at a given time
    of -- (instance) phase space variable
    wrt -- (instance) phase space variable, integration var

    outputs:
    F -- (ndarray, ndim = 1 or float) integrated result
    """

    return np.sum(f, axis = 1)*wrt.width

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
