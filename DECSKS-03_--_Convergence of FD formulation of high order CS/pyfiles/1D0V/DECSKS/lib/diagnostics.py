import DECSKS
import numpy as np
import numpy.linalg as LA
import sys # to retrieve smallest float for lower bound tolerance

def HighPrecisionE(number):
    """Converts a number into a string object
    while retaining a chosen degree of precision. This
    is designed to evade the truncation that is involved
    with str() so that outputs can store numbers with high
    precision

    inputs:
    number -- (number)

    outputs:
    string object with chosen precision in scientific notation
    """

    return "%.22e" % number

def calcs_and_writeout(sim_params,f,n,x,t):
    """orchestrates the calculation of various quantities, e.g.
    Lp norms, energy, electrostatic energy, ...

    inputs:
    files -- (dict) contains output filenames to be written to
    f -- (ndarray, ndim=3), f(t,x,v)
    n -- (int) time step number, t^n
    x -- (instance) space variable
    v -- (instance) velocity variable

    outputs:
    None
    """
    #I = "invariant", I1 = L1 norm invariant, etc.
    if sim_params['record_outputs'] == 'yes':
        I1 = L1(f,n,x)
        I2 = L2(f,n,x,t)
        S = entropy(f,n,x)
        # write to files

        sim_params['outfiles']['I1'].write(HighPrecisionE(I1) + '\n')
        sim_params['outfiles']['I2'].write(HighPrecisionE(I2) + '\n')
        sim_params['outfiles']['S'].write(HighPrecisionE(S) + '\n')

    if n == sim_params['Nt']:
        close_all_outfiles(sim_params)

    return I2

def L1(f,n,x):
    """computes the L1 norm

    inputs:
    f -- (ndarray, ndim=3), f(t,x,v)
    n -- (int) time step number, t^n
    x -- (instance) space variable

    outputs:
    I1 -- (float) L1 norm
    """

    return np.sum(f[n,:]) * x.width

def L2(f,n,x,t):
    """computes the square of the L2 norm. Note, the intended
    purpose of this computation is to compare with its deviation
    from the value at time zero. To minimize compounded errors
    from redundant operations, a squareroot is not taken here
    and should be applied later if desired,
    e.g. np.sqrt( (L2[t] - L2[0]) / L2[0])

    inputs:
    f -- (ndarray, ndim=3), f(t,x,v)
    n -- (int) time step number, t^n
    x -- (instance) space variable

    outputs:
    I2 -- (float) L2 norm
    """

    # compute the square of the L2 norm below to minimize
    # compounded error from repeated operations like squareroot
    # TODO this has been modified for convergence analysis temporarily, computes L2 of L1 error norm

    a, c = 3/4., 1/2.
    w_a, w_b, w_c = 0.03, 0.06, 0.1
    x_a, x_c = 0.25, -0.25

    f_exact = a * np.exp(-(((x.cells - 0.1) + x_a) / w_a) ** 2) + np.exp(-(((x.cells - 0.1)/w_b)**2)) + c*np.exp(-(((x.cells - 0.1) + x_c)/w_c) ** 2)

    error = f[n,:] - f_exact
    return np.sqrt(x.width / x.L) * LA.norm(error, 2)

def entropy(f,n,x):
    """computes the entropy S at time t^n,
        S = sum_i sum_j f_[n,i,j] * ln (f[n,i,j] + eps) dxdv

    inputs:
    f -- (ndarray, ndim=3), f(t,x,v)
    n -- (int) time step number, t^n
    x -- (instance) space variable
    v -- (instance) velocity variable

    outputs:
    S -- (float) entropy at time t^n
    """
    eps = sys.float_info.min # to evade taking np.log(0)
    return np.sum(f[n,:] * np.log(f[n,:] + eps)) * x.width

def close_all_outfiles(sim_params):
    """Closes all opened output files inside dictionary
    sim_params['outfiles']

    inputs:
    sim_params -- (dict) simulation parameters, includes dict of outfiles

    outputs:
    None
    """
    if sim_params['outfiles'] is not None:
        for outfile in sim_params['outfiles'].itervalues():
            outfile.close()
    return None
