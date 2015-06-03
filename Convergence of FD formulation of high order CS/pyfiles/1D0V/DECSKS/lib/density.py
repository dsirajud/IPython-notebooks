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
        f = np.zeros([t.N + 1, z1.N])
        f[0,:] = initial_profile(f[0,:], sim_params, z1)

    elif z2 is not None:
        f = np.zeros([t.N + 1, z1.N, z2.N])
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
        for i in range(x.N):
            for j in range(v.N):
                f0[i,j] =1 / np.sqrt(2*np.pi) * (1 + eps*np.cos(k*x.cells[i])) * np.exp(-v.cells[j] ** 2 / 2.)
        return f0

    if density == 'bump on tail':
        x,v = z1, z2
        for i in range(x.N):
            for j in range(v.N):
                f0[i,j] = 1 / np.sqrt(2*np.pi) * (1 + 0.04*np.cos(0.3*x.cells[i])) * ( 0.9*np.exp(-v.cells[j]**2 / 2.) + 0.2*np.exp(-4 * (v.cells[j] - 4.5) ** 2) )
        return f0

    if density == 'gaussian':
        mu, sigma = 0.0, 0.04
        return np.exp(-(z1.cells - mu)**2/(2*sigma**2))

    if density == 'n cosine bell':
        # 6th degree cosine bell, cos(2pix)**6 in |x| < 0.5
        n = 6
        for i in range(z1.N):
            if np.abs(z1.cells[i]) < 0.25:
                f0[i] = 0.1 + (np.cos(2*np.pi*z1.cells[i])) ** n
            else:
                f0[i] = 0.1
        return f0

    if density == 'rectangle and gaussian bell':

        for i in range(z1.N):
            if -0.4 <= z1.cells[i] <= -0.2:
                f0[i] = 0.1 + 1.0
            elif -0.1 <= z1.cells[i] <= 0.5:
                f0[i] = 0.1 + np.exp(- ((z1.cells[i] - 0.2)/0.04)**2)
            else:
                f0[i] = 0.1 + 0.0
        return f0

    if density == 'triangle':

        for i in range(z1.N):
            if np.abs(z1.cells[i]) < 0.25:
                f0[i] = 0.1 + (1 - 4*np.abs(z1.cells[i]))
            else:
                f0[i] = 0.1
        return f0

    if density == 'triple gaussian bell':

        f0[:] = 0.5*np.exp(-((z1.cells + 0.2) / 0.03)**2) + np.exp(-((z1.cells) / 0.06)**2) + 0.5*np.exp(-((z1.cells - 0.2) / 0.03)**2)
        return f0

    if density == 'quint gaussian bell':

        f0[:] = 0.5*np.exp(-((z1.cells + 0.25) / 0.03)**2) + 0.5*np.exp(-((z1.cells + 0.2) / 0.03)**2) + np.exp(-((z1.cells) / 0.06)**2) + 0.5*np.exp(-((z1.cells - 0.2) / 0.03)**2) + 0.5*np.exp(-((z1.cells - 0.25) / 0.03)**2)
        return f0

    if density == 'triple gaussian bell asymmetric':
        a, c = 3/4., 1/2.
        w_a, w_b, w_c = 0.03, 0.06, 0.1
        x_a, x_c = 0.25, -0.25
        f0[:]  = a * np.exp(-((z1.cells + x_a) / w_a) ** 2) + np.exp(-((z1.cells/w_b)**2)) + c*np.exp(-((z1.cells + x_c)/w_c) ** 2)
        return f0


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
    """
    ax, bx = sim_params['ax'], sim_params['bx']
    ne = single_integration(f[0,:,:], of = x, wrt = v)
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
