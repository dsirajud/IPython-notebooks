import DECSKS
import numpy as np
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

def calcs_and_writeout(sim_params,f,n,x,v):
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
        I1 = L1(f,n,x,v)
        I2 = L2(f,n,x,v)

        # electrostatic terms
        E = DECSKS.lib.fieldsolvers.Poisson(sim_params['ni'], f, x, v, n)
        IW = total_energy(f,n,x,v,E)
        WE = electrostatic_energy(x,E)
        S = entropy(f,n,x,v)
        # write to files

        sim_params['outfiles']['I1'].write(HighPrecisionE(I1) + '\n')
        sim_params['outfiles']['I2'].write(HighPrecisionE(I2) + '\n')
        sim_params['outfiles']['IW'].write(HighPrecisionE(IW) + '\n')
        sim_params['outfiles']['WE'].write(HighPrecisionE(WE) + '\n')
        sim_params['outfiles']['S'].write(HighPrecisionE(S) + '\n')

    if n == sim_params['Nt']:
        close_all_outfiles(sim_params)

    return None

def L1(f,n,x,v):
    """computes the L1 norm

    inputs:
    f -- (ndarray, ndim=3), f(t,x,v)
    n -- (int) time step number, t^n
    x -- (instance) space variable
    v -- (instance) velocity variable

    outputs:
    I1 -- (float) L1 norm
    """

    return np.sum(f[n,:,:]) * x.width * v.width

def L2(f,n,x,v):
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
    v -- (instance) velocity variable

    outputs:
    I2 -- (float) L2 norm
    """

    # compute the square of the L2 norm below to minimize
    # compounded error from repeated operations like squareroot
    return np.sum(f[n,:,:]**2) * x.width * v.width

def total_energy(f,n,x,v,E):
    """computes the total energy for a Vlasov-Poisson system
        IW = 1/2 sum_i sum_j f[n,i,j] dx dv + 1/2 sum_i E[i] dx

    inputs:
    f -- (ndarray, ndim=3), f(t,x,v)
    n -- (int) time step number, t^n
    x -- (instance) space variable
    v -- (instance) velocity variable
    E -- (ndarray, ndim=1), E(x) at t^n

    outputs:
    IW -- (float) total energy at time t^n in system
    """

    return 1/2.*np.sum(f[n,:,:] * v.cells **2) * x.width * v.width \
      + 1/2. * np.sum(E**2) * x.width

def electrostatic_energy(x,E):
    """computes the electrostic energy WE = 1/2 sum_i E[i] dx

    inputs:
    E -- (ndarray, ndim=1) E(x) at t^n
    x -- (instance) space variable

    outputs:
    WE -- (float) electrostatic energy at time t^n
    """
    return 1/2.* np.sum(E**2)* x.width

def entropy(f,n,x,v):
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
    return np.sum(f[n,:,:] * np.log(f[n,:,:] + eps)) * x.width * v.width

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
