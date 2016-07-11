import numpy as np
import numpy.linalg as LA
import DECSKS

#==============================================================================#
# ORCHESTRATORS
#==============================================================================#

def compute_electric_field_fd(fe, fi, x, vx, n, sim_params):

    phi = eval(sim_params['compute_electric_potential_phi_handle'][x.str])(fe, fi, x, vx, n, sim_params)

    # currently the finite difference weight matrix W_dn1 is a 6th order LTE to
    # match the 6th order LTE on the Poisson solve
    dphi = 1 / x.width * sim_params['W_dn1_LTE6'].dot(phi)
    Ex = -dphi

    return Ex

def compute_electric_field_fourier(fe, fi, x, vx, n, sim_params):
    """Computes self-consistent electric field E by solving Gauss' law
    using FFT/IFFT for two species.

    inputs:
    ni -- (float) uniform background density of ions,
                  in the future can take an input fi, to compute ni
    fe -- (ndarray, dim=2) electron density fe(x,v,n) at time step t^n
          shape = (x.Ngridpoints, vx.Ngridpoints)
    fi -- (ndarray, dim=2) electron density fi(x,v,n) at time step t^n
          shape = (x.Ngridpoints, vx.Ngridpoints)
    x -- (instance) spatial variable
    vx -- (instance) velocity variable
    n -- (int) time step number, t^n
    sim_params -- (dict) simulation parameters dictionary
        sim_params['xi'] -- (ndarray, ndim=1) wave number vector over x

    outputs:
    E -- (ndarray,dim=2) electric field, E(x,v) = E(x) at time t^n for all (i,j)
    """

    # here we solve for the field at each i = 0, 1, ..., x.Ngridpoints - 1

    # integrate int dv_x (fi(x,v_x) - fe(x,v_x)) = sum (fi[i,j] - fe[i,j]) * vx.width
    # where we over j = 0, 1, 2, ... , vx.Ngridpoints - 1. This is a midpoint
    # rule given our cell-centered grid
    n_total = np.sum(fi[:x.Ngridpoints-1, :vx.Ngridpoints-1] - fe[:x.Ngridpoints-1, :vx.Ngridpoints-1], axis = 1) * vx.width
    print "n_total.shape = "
    print n_total.shape

    # Fourier transform the total density, n_total, and the field E
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

#==============================================================================#
# SOLVERS
#==============================================================================#

#==============================================================================3
# TWO SPECIES SOLVERS
#==============================================================================#
# TODO: all of these only differ in their D matrix (already stored) and
# auxillary boundary value vector, phi_BC. We can group all these into a
# class of solvers rather than make a separate routine for all
#==============================================================================#

def Poisson_6th_PBC(fe, fi,
                x, vx, n,
                sim_params):
    """6th order LTE finite difference Poisson solver for periodic BCs

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
    # midpoint rule, integrate over vx space
    n_total = np.sum(fe - fi, axis = 1) * vx.width

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = B*n + phi_BC, phi_BC = 0

    # label the RHS as b = dx ** 2 * B*n
    b = x.width ** 2 * sim_params['Poisson_6th_order_FD_solver_matrices']['B'].dot(n_total)

    # solve D*phi = b
    phi = LA.solve(sim_params['Poisson_6th_order_FD_solver_matrices']['D'], b)

    # boundary conditions: periodic

    # periodic do not produce unique solutions but a family of solutions with arbitrary integration constant
    # that corresponds to a DC offset phi_avg, recenter so that phi_avg = 0

    phi_avg = np.sum(phi) * x.width / x.L
    phi -= phi_avg

    # generate the 2D map for every [i,j], note that each row is constant
    phi = np.outer(phi, np.ones([1,vx.N]))

    return phi

def Poisson_6th_LDBC_UDBC(fe, fi,
                x, vx, n,
                sim_params):
    """6th order LTE finite difference Poisson solver for

        - lower boundary: Dirichlet
        - upper boundary: Dirichlet

    see notebook DECSKS-04 for details on construction

    https://github.com/dsirajud/IPython-notebooks/blob/master/
    DECSKS-04%20--%20Design%20of%20a%206th%20order%20FD%20Poisson%20solver/
    DECSKS-04%20--%20Design%20of%20a%206th%20order%20FD%20Poisson%20solver.ipynb

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

    # midpoint rule, integrate over vx space
    n_total = np.sum(fe - fi, axis = 1) * vx.width

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = b
    #
    # where     b = x.width ** 2 * B*n + phi_BC

    # BOUNDARY CONDITIONS: already set up in lib.read in vector sim_params['phi_BC'][x.str]

    # RHS of Poisson equation
    b = x.width ** 2 * sim_params['Poisson_6th_order_FD_solver_matrices']['B'].dot(n_total) + sim_params['phi_BC'][x.str]

    # solve D*phi = b
    phi = LA.solve(sim_params['Poisson_6th_order_FD_solver_matrices']['D'], b)

    # generate the 2D map for every [i,j], note that each row is constant
    phi = np.outer(phi, np.ones([1,vx.N]))

    return phi

def Poisson_6th_LNBC_UDBC(fe, fi,
                x, vx, n,
                sim_params):
    """6th order LTE finite difference Poisson solver for

        - lower boundary: Neumann
        - upper boundary: Dirichlet

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
    # midpoint rule, integrate over vx space
    n_total = np.sum(fe - fi, axis = 1) * vx.width

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = b
    #
    # where     b = x.width ** 2 * B*n + phi_BC

    # BOUNDARY CONDITIONS:

    # the upper DBC has been read in and stored in sim_params['phi_BC'][x.str][-1] based on params_boundaryconditions.dat input during lib.read execution
    # the lower NBC can be either a symmetric condition (boundary value on derivative = 0) or self-consistent (boundary value depends on wall charge)
    if sim_params['BC']['phi'][x.str]['lower'] == 'SYMMETRIC' or sim_params['BC']['phi'][x.str]['lower'] == 'SYMMETRY':
        pass # phi_BC[x.str][0] = 0 , which is the symmetry BC on phi is already zero from lib.read initialization of the ndarray
    elif sim_params['BC']['phi'][x.str]['lower'] == 'SELF-CONSISTENT':
        # TODO check sign on this, probably dphi(a) = -x.width * sigma, and dphi(b) = +x.width * sigma
        sim_params['phi_BC'][x.str][0] = -x.width * sim_params['sigma'][x.str]['lower']

    b = x.width ** 2 * sim_params['Poisson_6th_order_FD_solver_matrices']['B'].dot(n_total) + sim_params['phi_BC'][x.str]

    # solve D*phi = b
    phi = LA.solve(sim_params['Poisson_6th_order_FD_solver_matrices']['D'], b)

    # generate the 2D map for every [i,j], note that each row is constant
    phi = np.outer(phi, np.ones([1,vx.N]))

    return phi


def Poisson_6th_LDBC_UNBC(fe, fi,
                x, vx, n,
                sim_params):
    """6th order LTE finite difference Poisson solver for

        - lower boundary: Dirichlet
        - upper boundary: Neumann

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

    # not implemented yet, a NotImplementedYet error would have already been thrown in lib.read

    return None

def Poisson_6th_LDBC_LNBC(fe, fi,
                x, vx, n,
                sim_params):
    """6th order LTE finite difference Poisson solver for LDBC/LNBC

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
    # not implemented yet, a NotImplementedYet error would have already been thrown in lib.read

    return None

def Poisson_6th_UDBC_UNBC(fe, fi,
                x, vx, n,
                sim_params):
    """6th order LTE finite difference Poisson solver for LDBC/LNBC

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
    # not implemented yet, a NotImplementedYet error would have already been thrown in lib.read

    return None

#------------------------------------------------------------------------------#
# ONE SPECIES SOLVERS (one evolved species + one const. background)
# 1S = one species
#------------------------------------------------------------------------------#
def Gauss_1S(ni, f, x, vx, n, sim_params):
    """Computes self-consistent electric field E by solving Gauss' law
    using FFT/IFFT.

    ** This routine calculates the E field and copies it over
       vx.N columns so that each [i,j] gives the field directly at
       such a point so that the velocities at every [i,j] can be
       updated straightforwardly.

       This routine could be used to calculate electric field related quantities
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

def Poisson_6th_PBC_1S(ni, f,
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

#==============================================================================#
# MISCELLANEOUS METHODS USED ABOVE
#==============================================================================#
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
    if we wish to integrate along rows we specifiy axis = 0
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
