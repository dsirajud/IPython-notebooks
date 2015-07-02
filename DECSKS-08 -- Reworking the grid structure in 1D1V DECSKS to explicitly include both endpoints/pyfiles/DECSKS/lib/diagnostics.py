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
        I2 = L2(f,n,x)
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

def L2(f,n,x):
    """computes the L2 norm of the error.

    inputs:
    f -- (ndarray, ndim=3), f(t,x,v)
    n -- (int) time step number, t^n
    x -- (instance) space variable

    outputs:
    I2 -- (float) L2 norm
    """
    f_exact = f[0,:]
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
