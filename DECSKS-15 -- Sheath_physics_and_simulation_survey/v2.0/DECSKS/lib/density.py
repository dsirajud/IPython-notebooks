import numpy as np

def setup(sim_params, t, z1, z2 = None):
    """Returns an initial density

    inputs:
    sim_params -- (dict) simulation parameters
    t -- (instance) time
    z1 -- (instance) phase space variable, e.g. x
    z2 -- (instance) phase space variable, e.g. y, v

    outputs:
    f -- (ndarray, ndim = 2,3) density f with f(x,v = None, t)
         loaded
    """

    if z2 is None:
        # init density container with dims z1.N and zero time + t.N steps
        f = np.zeros([t.Ngridpoints, z1.Ngridpoints])
        f[:,0] = initial_profile(f[0,:], sim_params, z1)

    elif z2 is not None:
        f = np.zeros([t.Ngridpoints, z1.Ngridpoints, z2.Ngridpoints])
        f[0,:,:] = initial_profile(f[0,:,:], sim_params, z1, z2)

    return f

def initial_profile(f0, sim_params, z1, z2 = None):
    """Returns the initial density specified in input file

    inputs:
    f0 -- (ndarray, ndim = 1,2) density container at t = 0
    density -- (str)
    sim_params -- (dict) simulation parameters
    z1 -- (instance) phase space variable
    z2 -- (instance) phase space variable

    outputs:
    f0 -- (ndarray, ndim = 1,2) initial density according
          to sim_params['density']
    """
    density = sim_params['density']
    if density == 'landau':
        x,v = z1, z2
        eps, k = 0.01, 0.5
        print "initializing Landau profile"
        for i in range(x.Ngridpoints):
            for j in range(v.Ngridpoints):
                f0[i,j] =1 / np.sqrt(2*np.pi) * (1 + eps*np.cos(k*x.gridvalues[i])) * np.exp(-v.gridvalues[j] ** 2 / 2.)
        return f0

    if density == 'bump on tail':
        x,v = z1, z2
        for i in range(x.N):
            for j in range(v.N):
                f0[i,j] = 1 / np.sqrt(2*np.pi) * (1 + 0.04*np.cos(0.3*x.gridvalues[i])) * ( 0.9*np.exp(-v.gridvalues[j]**2 / 2.) + 0.2*np.exp(-4 * (v.gridvalues[j] - 4.5) ** 2) )
        return f0

    if density == 'gaussian':
        mu, sigma = 0.0, 0.04
        f0[:] = np.exp(-(z1.gridvalues - mu)**2/(2*sigma**2))
        return f0

    if density == 'n cosine bell':
        # 6th degree cosine bell, cos(2pix)**6 in |x| < 0.5
        n = 6
        for i in range(z1.N):
            if np.abs(z1.gridvalues[i]) < 0.25:
                f0[i] = 0.1 + (np.cos(2*np.pi*z1.gridvalues[i])) ** n
            else:
                f0[i] = 0.1
        return f0

    if density == 'rectangle and gaussian bell':

        for i in range(z1.N):
            if -0.4 <= z1.gridvalues[i] <= -0.2:
                f0[i] = 0.1 + 1.0
            elif -0.1 <= z1.gridvalues[i] <= 0.5:
                f0[i] = 0.1 + np.exp(- ((z1.gridvalues[i] - 0.2)/0.04)**2)
            else:
                f0[i] = 0.1 + 0.0
        return f0

    if density == 'triangle':

        for i in range(z1.N):
            if np.abs(z1.gridvalues[i]) < 0.25:
                f0[i] = 0.1 + (1 - 4*np.abs(z1.gridvalues[i]))
            else:
                f0[i] = 0.1
        return f0

    if density == 'triple gaussian bell':

        f0[:] = 0.5*np.exp(-((z1.gridvalues + 0.2) / 0.03)**2) + np.exp(-((z1.gridvalues) / 0.06)**2) + 0.5*np.exp(-((z1.gridvalues - 0.2) / 0.03)**2)
        return f0

def global_conservation_check(sim_params, f_new, n):
    """Checks if mass is conserved over the remapping procedure

    NOTE: mass conservation is intrinsic to the algorithm. If you
    see mass being lost or added, the implemenetation was not done
    correctly.

    inputs:
    f_new -- (ndarray, ndim = 2) container with remapped MC density
    f_old -- (ndarray, ndim = 2) density from previous time step
    n -- (int) current time step that f_new pertains to, used
         to alert user to which time step any non-conservation occurred

    outputs:
    None
    """

    TOL = 1.0e-14

    # the key/value pair sim_params['m_0'] is assigned in main.py after
    # the density profile f is instantiated
    mass_difference = np.abs( sim_params['m_0'] - np.sum(np.abs(f_new)))
    if mass_difference > TOL:
        print 'mass difference = %.40e, density is not conserved \
        globally at time step = %d' % (mass_difference, n)
    return None

def conservation_check(f_new, f_old, i, n):
    """Checks if mass is conserved from the remap step from one
    MC at f_old[i] to contiguous cells in f_new container

    inputs:
    f_new -- (ndarray, ndim = 1) container with remapped MC density
    f_old -- (ndarray, ndim = 1) density from previous time step
    i -- (int) prepoint of MC
    n -- (int) current time step

    outputs:
    None
    """
    eps = 1.0e-14
    mass_difference = np.abs( sum(f_new) - f_old[i])
    if mass_difference > eps:
        print 'mass difference = %.40e, density is not conserved \
        locally at time step = %d, from prepoint %i' % (mass_difference, n, i)
    return None

def cold_background(f,x,v,sim_params):
    """Returns a cold background density according to initial
    density f[x,v,0] per quasineutrality

    inputs:
    f -- (ndarray, ndim = 3) f(x,v,t) for all times
    x -- (instance) spatial variable
    v -- (instance) velocity variable
    sim_params -- (dict) simulation parameters

    outputs:
    ni -- (float) constant ion background density

    Note: for periodic BCs, z.N != z.Ngridpoints, hence
    we pass only the active gridpoints [0:z.N].
    for all other BCs, z.N = z.Ngridpoints  so the following
    pass is general
    """
    ax, bx = sim_params['ax'], sim_params['bx']
    ne = single_integration(f[0,:x.N,:v.N], of = x, wrt = v)
    Ne = single_integration(ne, wrt = x)
    ni = Ne / (bx - ax)
    return ni

def single_integration(f, of = None, wrt = None):
    """integrates once a single variable or two variable
    function, i.e. computes integral f(z,wrt) d(wrt) = f(z),
    or integral f(wrt) d(wrt) = F. For the case of a two
    variable function, 'of' is the unintegrated variable
    such that f is a function 'of' that variable after it was
    integrated with respect to the variable 'wrt'.

    inputs:
    f -- (ndarray, ndim = 1,2) density at a given time
    of -- (instance) phase space variable
    wrt -- (instance) phase space variable, integration var

    outputs:
    F -- (ndarray, ndim = 1 or float) integrated result
    """
    z = of
    if z is not None:
        F = np.zeros(z.N)
        for i in z.prepoints:
            F[i] = riemann_sum(f[i,:],wrt)
    elif z is None:
        F = sum(f)*wrt.width
    return F

def riemann_sum(f, wrt):
    """Computes integral f(wrt) d(wrt) when spacing wrt.width
    is uniform on mesh

    inputs:
    f -- (ndarray, ndim = 1) 1D array
    wrt -- (instance) phase space variable

    outputs:
    ne -- (float) integrated result
     """
    ne = sum(f)*wrt.width
    return ne

    # if non-uniform spacing:
    #z = wrt
    #ne = 0
    #for i in range(z.N):
    #    ne += f[i]*z.width
