import numpy as np
from math import factorial

def fourier(d, f, z, sim_params):
    wave_index = np.arange(z.N)
    xi = np.where(wave_index <= z.N / 2,
                  2*np.pi*wave_index / z.L,
                  2*np.pi*(wave_index - z.N) / z.L)

    # for x advection
    if z.str == 'x':
        # generate a 2D matrix appropriate for numpy multiplication with the 2D array Ff below
        xi = np.outer(xi, np.ones(eval(sim_params['phasespace_vars'][1]).N))

    elif z.str[0] == 'v':
        # generate a 2D matrix appropriate for numpy multiplication with the 2D array Ff below
        xi = np.outer(xi, np.ones(eval(sim_params['phasespace_vars'][0]).N))

    for q in range(1,sim_params['N']):
        d[q,:,:] = trigonometric2D(f, z, q, xi)

def trigonometric2D(f, z, q, xi):
    """Computes derivatives of density for derivative coeffs d in FN methods

    inputs:
    f -- (ndarray, dim=1) f(z1,z2=const,z3=const,..., t = n-1)
    z -- (instance) phase space variable being convected
    q -- (int) order of desired derivative
    xi -- (ndarray, dim=2) 1D wave numbers stacked in columns for column-wise
          fourier transform.
    sim_params -- (dict) simulation parameters
    K -- (ndarray, ndim=1) Windowed Fourier transform kernel

    outputs:
    d -- (ndarray, dim=1) qth derivative coefficient for all z.prepoints
    """

    Ff = np.fft.fft(f)

    # zero out negligible contributions
    A = max(Ff)
    Ff_min = A*(2.0e-15)
    Ff = np.where(Ff < Ff_min, 0, Ff)

    D = (1j*xi*z.width) ** q * Ff
    d = np.real(np.fft.ifft(D))

    return d

def trigonometric(f, z, q, sim_params, K = None):
    """Computes derivatives of density for derivative coeffs d in FN methods

    inputs:
    f -- (ndarray, dim=1) f(z1,z2=const,z3=const,..., t = n-1)
    z -- (instance) phase space variable being convected
    q -- (int) order of desired derivative
    sim_params -- (dict) simulation parameters
    K -- (ndarray, ndim=1) Windowed Fourier transform kernel

    outputs:
    d -- (ndarray, dim=1) qth derivative coefficient for all z.prepoints
    """
    xi = np.zeros(z.N)
    for r in range(z.N):
        if r <= z.N/2 :
            xi[r] = 2*np.pi*r / z.L
        else:
            xi[r] = 2*np.pi*(r - z.N) / z.L

    Ff = np.fft.fft(f)

    if sim_params['WindowedFilter'] == 'YES':
        FK = np.fft.fft(K) # used to be np.real(np.fft.fft(K))
        Ff = FK*Ff

    A = max(Ff)
    Ff_min = A*(2.0e-15)
    for r in range(z.N):
        if np.abs(Ff[r]) < Ff_min:
            Ff[r] = 0

    D = (1j*xi*z.width) ** q * Ff
    d = np.real(np.fft.ifft(D))

    return d

def finite_differences_6th_order_periodic_BCs(f,z):
    """Computes the derivative coefficients for FD methods in high order CS
    for periodic BCs (hence the use of np.mod())

    inputs:
    f -- (ndarray, ndim=1) f(z1,z2=const,z3=const,..., t = n-1)
    z -- (instance) phase space variable being convected

    outputs:
    d -- (ndarray, ndim=2) d[:,q] is the qth derivative for all points i
        weighted by a power (z.width ** q) per CS correction prescription
    """
    d = np.zeros([z.N,5])

    for i in z.prepoints:
        d[i,1] = (1/12.)*(f[np.mod(i-2,z.N)] - 8*f[np.mod(i-1,z.N)] + 8*f[np.mod(i+1,z.N)] - f[np.mod(i+2,z.N)])
        if 50 < i < 100:
            print "d[%d,1] = %g" % (i, d[i,1])
        d[i,2] = (1/12.)*(-f[np.mod(i-2,z.N)] + 16*f[np.mod(i-1,z.N)] - 30*f[np.mod(i,z.N)] + 16*f[np.mod(i+1,z.N)] - f[np.mod(i+2,z.N)])
        d[i,3] = (1/2.)*(-f[np.mod(i-2,z.N)] + 2*f[np.mod(i-1,z.N)] - 2*f[np.mod(i+1,z.N)] + f[np.mod(i+2,z.N)])
        d[i,4] = f[np.mod(i-2,z.N)] - 4*f[np.mod(i-1,z.N)] + 6*f[np.mod(i,z.N)] - 4*f[np.mod(i+1,z.N)] + f[np.mod(i+2,z.N)]
    return d

def fd(d, f, z, sim_params):
    W = sim_params['W']
    Wz = W[z.str]

    for dn in range(1,sim_params['N']):
        d[dn,:,:] = Wz[dn,:,:].dot(f)
    return d

def finite_differences(f, z, sim_params):
    """Computes derivatives via explicit finite differences for each grid point
    i in a bounded domain. FD schemes here have LTE = O(dz^6) which permit
    correcting CS up to order 5 (global error) using the derivatives
    dn = 1, 2, 3, 4. The philosophy followed here is for derivatives to be
    calculated "as centered as reasonably achievable" (ACARA), hence
    near edges the schemes shift their stencils only when necessary
    and do so gridpoint-by-gridpoint in order to capture information from
    as many grid points on either side as possible

    inputs:
    f -- (ndarray, ndim=1) f(z1,z2=const,z3=const,..., t = n-1)
    z -- (instance) phase space variable being convected
    sim_params -- (dict) simulation parameters (contained FD_schemes dictionary)

    outputs:
    d -- (ndarray, ndim=2) d[:,q] is the qth derivative for all points i
        weighted by a power (z.width ** q) per CS correction prescription
    """
    #p = sim_params['N'] + 1
    d = np.zeros([z.N, sim_params['N']])
    imax = z.N - 1

    for dn in range(1,sim_params['N']): # dn = "derivative number", dn-th derivative

        # retrieve subdictionary of all FD schemes for derivative number dn
        FD_schemes = sim_params['FD_schemes']['dn' + str(dn)]
        p = sim_params['N'] - dn + 1 # = GE + 1 = LTE, truncation error
        # i.e. the requirement of p decreases with each dn if using the tables
        # generated by ./bin/finite_difference_schemes/
        # /generate_tables_of_finite_difference_schemes_required_for_a_chosen_GE_on_CS.py

        stencil_size = p + dn
        stencil_center = stencil_size // 2

        for i in z.prepoints:
            if i < stencil_center: # left-edge of domain
                d[i,dn] = FD(f, z, FD_schemes,
                             handedness = 'forward',
                             dn = dn,
                             z0 = i,
                             asymmetry = i
                             )

            elif (imax - i) < stencil_center: # right-edge of domain
                d[i,dn] = FD(f, z, FD_schemes,
                             handedness = 'backward',
                             dn = dn,
                             z0 = i,
                             asymmetry = imax - i
                             )

            elif i >= stencil_center and (imax - i) >= stencil_center: # interior of domain
                if np.mod(stencil_size, 2) == 1: # stencil size is odd
                    #                    print "stencil size is odd, using central scheme, dn = %d" % dn
                    d[i,dn] = FD(f, z, FD_schemes,
                                 handedness = 'central',
                                 dn = dn,
                                 z0 = i,
                                 asymmetry = 0)

                else: # stencil is even (no center, and no 'central' scheme)
                    #                    print "stencil size is even, using forward - %d, dn = %d" % (stencil_center - 1, dn)
                    # stencil_center - 1 defines the greatest asymmetry
                    # on either scheme (i.e. closest to central differencing)

                    d[i,dn] = FD(f, z, FD_schemes,
                                 handedness = 'forward', # or 'backward'
                                 dn = dn,
                                 z0 = i,
                                 asymmetry = stencil_center - 1
                                 )

    return d

def FD(f, z, FD_schemes,
       handedness = None,
       dn = None,
       z0 = None,
       asymmetry = None
       ):
    """Computes derivatives via finite differences for a grid point z0
    in a bounded domain by selecting one scheme among the family of
    explicit schemes of order LTE = O(z^{dn}). Schemes are labeled according to:

    'forward'

        asymmetry:

        0 = forward differencing, all grid points used are forward of i
        1 = 1 grid points is backward of i, rest are forward
        2 = 2 grid points are backward of i, rest are forward,
        ...

    'central'

        asymmetry:

        0 = central differencing, symmetric sampling of grid points about i

    'backward'

        asymmetry:

        2 = 2 grid points forward of i, rest are backward
        1 = 1 grid point forward of i, rest are backward
        0 = backward differencing,
        ...

    difference weights are calculated from
    ./bin/finite_difference_schemes/
    generate_finite_difference_schemes_required_for_a_chosen_GE_on_CS.py
    routine and stored in etc/finite_difference_schemes/ as

        Table_of_difference_schemes.dat
        f1_FD_coefficients.dat <-- first derivative schemes
        f2_FD_coefficients.dat <-- second derivative schemes
        ...

    the single derivative files are what are read in and stored in a dictionary
    FD_schemes (sim_params['FD_schemes']) in lib/read.py

    FD_schemes -- (dict) each key permits selection of the derivative number
                  'dn#' whose values are themselves subdictionaries that have
                  keys indicating the handedness of the scheme: (['forward'],
                  ['backward'], ['central']) whose values are themselves
                  subdictionaries that have keys corresponding to the degree of
                  asymmetry in in the stencil sampling '0', '1', '2', ..., and
                  each of these have values that are themselves subdictionaries
                  whose keys correspond to values that are lists 'w' and
                  'stencil' which give the difference weights and stencil for
                  the selected scheme

    for example, to retrieve the dictionary characterizing a forward scheme
    with a degree of asymmetry = 2 on the third derivative,
    we extract the required dictionary first

    # dictionary of all schemes
    FD_schemes = sim_params['FD_schemes']

    # dictionary for specific scheme
    FD_scheme = FD_schemes['dn' + str(3)]['forward']['2']

    # then, store the specific weights w and stencil used according to
    w = FD_scheme['w']
    stencil = FD_scheme['stencil']

    inputs:
    f -- (ndarray, ndim=1) f(z1,z2=const,z3=const,..., t = n-1)
    z -- (instance) phase space variable being convected
    FD_schemes -- (dict) all FD schemes for derivative number dn
    handedness -- (str) 'forward','backward','central' of sampled grid pts
    dn -- (int) derivative number, e.g. dn = 2 = second derivative
    z0 -- (int) grid point "i" where derivative is to be calculated
    asymmetry -- (int) degree of asymmetry in sampling, e.g. 0, 1, 2

    outputs:
    d -- (float) derivative coefficient in CS algorithm (i.e. z.width ** dn * dnf)
         where dnf is the dn'th derivative of f at grid point i

         Note: calcualting d *= 1 / z.width **n after the loop would give the
         actual derivative, the d coefficient in high order
         CS is deriv * z.width ** n
    """
    # retrieve specific dictionary pertaining to specific scheme according
    # to handedness and degree of asymmetry
    FD_scheme = FD_schemes[handedness][str(asymmetry)]

    w = FD_scheme['w']
    stencil = FD_scheme['stencil']

    d = 0
    for s in range(len(stencil)):
        d += w[s]*f[z0 + stencil[s]]

    return d

def assemble_finite_difference_weight_matrices(sim_params, x, v):
    """Assembles the dictionary W that contains both finite diffference
    weight matrices for each phase space variable

    inputs:
    sim_params -- (dict) simulation parameters
    x -- (instance) configurational variable
    v -- (instance) velocity variable

    outputs:
    W -- (dict) W['x'] and W['v'] give the weight matrices
         W['x'][dn,x.N,x.N] and W['v'][dn,v.N,v.N] where each
         2D matrix W['z'][dn,:,:] gives the weight matrix for the
         dn-th derivative of the variable z
    """

    W = {}
    W['x'] = assemble_finite_difference_weight_matrix(sim_params, x)
    W['vx'] = assemble_finite_difference_weight_matrix(sim_params, v)

    return W

def assemble_finite_difference_weight_matrix(sim_params, z):
    """Assembles a matrix corresponding to the weights of in
    the finite difference computation of derivatives, i.e.
    assembles the weight matrix W, in

        1 / x.width * Wf =  df,     i.e. W are the difference
                                    coefficients, which do not
                                    contain the width of the
                                    abscissa value, e.g. x.width

    where f and df are vectors of length z.N in the 1D case.

    inputs:
    sim_params -- (dict) simulation parameters
    z -- (instance) phase space variable

    outputs:
    Wz -- (ndarray, ndim=3) Wz[dn, z.N, z.N] where each 2D matrix Wz[dn,:,:]
          is the weight matrix W corresponding to the dn-th derivative
    """
    imax = z.N - 1
    Wz = np.zeros([sim_params['N'], z.N, z.N])

    for dn in range(1,sim_params['N']):

        W_dn = np.zeros([z.N, z.N])
        p = sim_params['N'] - dn     # LTE of scheme on dn-th derivative so that
                                     #  LTE[z.width ** dn * dnf] = O(N)

        FD_schemes = sim_params['FD_schemes']['dn' + str(dn)]
        stencil_size = p + dn
        stencil_center = stencil_size // 2

        for i in range(z.N):

            if i < stencil_center:
                handedness = 'forward'
                asymmetry = str(i)
            elif imax - i < stencil_center:
                handedness = 'backward'
                asymmetry = str(imax - i)
            else:
                if np.mod(stencil_size,2) == 1:
                    handedness = 'central'
                    asymmetry = str(0)
                else:
                    handedness = 'forward'
                    asymmetry = str(stencil_center - 1)

            FD_scheme = FD_schemes[handedness][asymmetry]
            w = FD_scheme['w']
            stencil = FD_scheme['stencil']

            W_dn[i, i + np.array(stencil)] = w # load all weights at once into W_dn

        Wz[dn,:,:] = W_dn

    return Wz

def assemble_finite_difference_weight_matrix_single_derivative(sim_params, z, dn = 1, LTE = 6):
    """Assembles a matrix corresponding to the weights of in
    the finite difference computation of derivatives, i.e.
    assembles the weight matrix W, in

        Wf = df

    where f and df are vectors of length z.N.

    inputs:
    sim_params -- (dict) simulation parameters
    z -- (instance) phase space variable

    outputs:
    W -- (ndarray, ndim=3) Wz[dn, z.N, z.N] where each 2D matrix Wz[dn,:,:]
          is the weight matrix W corresponding to the dn-th derivative
    """
    imax = z.N - 1
    W = np.zeros([z.N, z.N])

    FD_scheme = sim_params['FD_scheme_dn1']['dn' + str(dn)]['LTE' + str(LTE)]
    stencil_size = LTE + dn
    stencil_center = stencil_size // 2

    for i in range(z.N):
        if i < stencil_center:
            handedness = 'forward'
            asymmetry = str(i)
        elif imax - i < stencil_center:
            handedness = 'backward'
            asymmetry = str(imax - i)
        else:
            if np.mod(stencil_size,2) == 1:
                handedness = 'central'
                asymmetry = str(0)
            else:
                handedness = 'forward'
                asymmetry = str(stencil_center - 1)

        w = FD_scheme[handedness][asymmetry]['w']
        stencil = FD_scheme[handedness][asymmetry]['stencil']

        W[i, i + np.array(stencil)] = w # load all weights at once into W_dn

    return W