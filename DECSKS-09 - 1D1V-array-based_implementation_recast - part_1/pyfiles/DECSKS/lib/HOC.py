import numpy as np
import scipy
import numpy.ma as ma

def Beta_matrix(sim_params,
                alpha,
                z):
    """constructs the B matrix, whose columns are
    the beta vectors (shape = N x 1) for each value
    of the generalized velocity for the advecting
    variable z.

    inputs:
    sim_params -- (dict) simulation parameters
    alpha -- (ndarray, ndim=2) CFL.frac matrix, shape = (z1.N, z2.N)"""

    # local copies of A matrices
    A_neg = sim_params['A_matrix']['-1']
    A_pos = sim_params['A_matrix']['1']

    N_arr = np.outer(np.arange(1,sim_params['N']+1), np.ones([1,sim_params['N' + z.v_str]]))

    # for broadcasting operations below to an N x v.N ("v" is the generalized velocity),
    # we require an identically dimensioned matrix. Hence, we slice the rows up to N values.
    # There is no loss of information here given that every row entry in a given column is constant.
    
    alpha_hat = alpha[:sim_params['N'],:sim_params['N' + z.v_str]]
    alpha_tilde = ma.array(alpha_hat ** N_arr / scipy.misc.factorial(N_arr))

    # mask negative values
    alpha_tilde[alpha_hat < 0] = ma.masked

    # operate on only positive values
    beta_neg = ma.dot(A_pos, alpha_tilde)

    # invert mask to mask all nonnegative values
    alpha_tilde.mask = np.logical_not(alpha_tilde.mask)

    # operate on only negative values
    beta_pos = ma.dot(A_neg, alpha_tilde)

    # consolidate all columns in a single matrix
    B = np.zeros([sim_params['N'], sim_params['N' + z.v_str]])

    # wherever beta_neg.mask is False (i.e. unmasked value), assign beta_neg, otherwise beta_pos
    B = np.where(beta_neg.mask == False, beta_neg.data, beta_pos.data)

    return B

def correctors(sim_params, CFL, z):
    """computes the correction coefficients c

    inputs:
    sim_params -- (dict) simulation parameters
    B -- (ndarray, ndim=2) matrix containing column
         vectors for each beta = beta(z_opp)
         in DECSKS-09 these are given symbols Bx and Bv

         shape = (N, z_adv.N)

    outputs:
    c -- (ndarray, ndim=2) correction coefficients
         with shape = (N, z_notadv.N)
         where z_notadv means the not advecting phase
         space variable in the scope of a 2D advection
         implementation
    """
    I_alternating = np.diag( (-np.ones(sim_params['N']))  ** np.arange(sim_params['N']) )

    B = Beta_matrix(sim_params,
                     CFL.frac,
                     z)

    c = I_alternating.dot(B)
    return c

def beta(p,a):
    """Returns the pth coefficient beta.

    inputs:
    p -- (int) pth coefficient
    a -- (float) uncorrected fractional distance of MC

    output:

    beta_p -- pth beta coefficient for given vel. or accel.
    """

    beta_p = a ** (p+1) / factorial(p+1)

    if p > 0:

        for q in range(p):

            if a >= 0:

                beta_p -= beta(q,a) / factorial(p + 1 - q)

            else: # a < 0

                beta_p -= (-1) ** (p+q) * beta(q,a) / factorial(p + 1 - q)

    return beta_p
# ........................................................................... #
# DEPRECATED in v2.0

def beta_m(a, B, N):
    A = np.zeros([N,N])
    alpha = np.zeros([N,1])
    for i in range(N):
        alpha[i,0] = a ** (i+1) / factorial(i+1)
        for j in range(i+1):
            A[i,j] = B[i-j]/factorial(i-j)
    return A.dot(alpha)
# ........................................................................... #
def kernel(z):
    """Returns a Gaussian enveloped window filtered wave number
    in Fourier space.

    inputs:
    z -- (instance) phase space variable

    output:
    K -- filtered wave number K(xi)
    """
    sigma = 4
    mu    = 0.0

    K = np.sinc(z.cells /z.width) * np.exp(- (z.cells - mu)**2 / (2*sigma**2))

    return K
# ........................................................................... #
