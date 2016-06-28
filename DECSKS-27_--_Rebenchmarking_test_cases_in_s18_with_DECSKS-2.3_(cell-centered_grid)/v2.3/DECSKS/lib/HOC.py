import numpy as np
import scipy
import numpy.ma as ma

def Beta_matrix(sim_params, z, vz):
    """constructs the B matrix, whose columns are
    the beta vectors (shape = N x 1) for each value
    of the generalized velocity for the advecting
    variable z.

    See DECSKS-09 part 1 for details of this calculation

    inputs:
    sim_params -- (dict) simulation parameters
    z.CFL.frac -- (ndarray, ndim=2) contains the fractional CFL numbers
                  for every [i,j]


    outputs:

    B -- (ndarray, ndim=2), shape = (N, vz.N)

            for x-advection: B.shape = (N, vx.N)
            for v-advection: B.shape = (N, x.N), note that vz.N = ax.N = x.N here
    """

    # local copies of A matrices
    A_neg = sim_params['A_matrix']['-1']
    A_pos = sim_params['A_matrix']['1']

    N_arr = np.outer(np.arange(1,sim_params['N']+1), np.ones([1,vz.N]))

    # for broadcasting operations below to an N x vz.N array
    # we require an identically dimensioned matrix. Hence, we slice the rows up to N values.
    # There is no loss of information here given that every row entry in the z.CFL.frac matrix
    # is constant, given z.CFL.frac = z.CFL.frac (vz.prepointvaluemesh)

    # Naming per DECSKS-09 notebook:
    #
    #        alpha = z.CFL.frac, shape: (z.N, vz.N)
    #        alpha_hat = truncated z.CFL.frac[:N, :vz.N], vz.N = z.CFL.frac.shape[1]
    #        alpha_tilde : alpha_tilde[q,j] = alpha_hat ** q,q = 0, 1, ... N-1
    #

    # catch cases where N > z.N, need to append extra (N - z.N)
    # copies of any row (all the same) to the remaining rows of alpha_hat
    if z.CFL.frac.shape[0] < sim_params['N']:
        alpha_hat = np.zeros([sim_params['N'], z.CFL.frac.shape[1]])
        alpha_hat[:z.CFL.frac.shape[0],:] = z.CFL.frac
        N_extra_rows = sim_params['N'] - z.CFL.frac.shape[0]

        # generate the object of the required size to fit the remaining submatrix
        # in alpha_hat not already filled, compute an outer product with any row (all same)
        # here, we take the first row for obviousness
        alpha_hat_extras = np.outer( np.ones(N_extra_rows), z.CFL.frac[0,:])
        alpha_hat[z.CFL.frac.shape[0]:z.CFL.frac.shape[0] + N_extra_rows,:] = alpha_hat_extras

    else:
        alpha_hat = z.CFL.frac[:sim_params['N'],:z.CFL.frac.shape[1]]

    alpha_tilde = ma.array(alpha_hat ** N_arr / scipy.misc.factorial(N_arr))

    mask_neg = (alpha_hat < 0)
    mask_pos = np.logical_not(mask_neg)

    # mask out negative values, leave only positives (>= 0)
    alpha_tilde.mask = mask_neg

    # operate on only positive values (>= 0)
    beta_pos = ma.dot(A_pos, alpha_tilde)

    # mask out positive (>= 0), leave only negatives
    alpha_tilde.mask = mask_pos

    # operate on only negative values
    beta_neg = ma.dot(A_neg, alpha_tilde)

    # consolidate all columns in a single matrix
    B = np.zeros([sim_params['N'], vz.N])

    # wherever beta_neg.mask is False (i.e. unmasked value), assign beta_neg, otherwise beta_pos
    B = np.where(mask_neg == True, beta_neg.data, beta_pos.data)

    return B

def Beta_matrix_for_given_stage(sim_params, z, vz, stage):
    """constructs the B matrix, whose columns are
    the beta vectors (shape = N x 1) for each value
    of the generalized velocity for the advecting
    variable z.

    See DECSKS-09 part 1 for details of this calculation

    inputs:
    sim_params -- (dict) simulation parameters
    z.CFL.frac -- (ndarray, ndim=2) contains the fractional CFL numbers
                  for every [i,j]


    outputs:

    B -- (ndarray, ndim=2), shape = (N, vz.N)

            for x-advection: B.shape = (N, vx.N)
            for v-advection: B.shape = (N, x.N), note that vz.N = ax.N = x.N here
    """

    # for broadcasting operations below to an N x vz.N array
    # we require an identically dimensioned matrix. Hence, we slice the rows up to N values.
    # There is no loss of information here given that every row entry in the z.CFL.frac matrix
    # is constant, given z.CFL.frac = z.CFL.frac (vz.prepointvaluemesh)

    # Naming per DECSKS-09 notebook:
    #
    #        alpha = z.CFL.frac, shape: (z.N, vz.N)
    #        alpha_hat = truncated (or extended) z.CFL.frac[:N, :vz.N] up to N rows
    #                    note that all rows are identical so copying extra rows
    #                    should N > z.N (hence, z.CFL.frac.shape[0] = z.N) is
    #                    performed at no consequence or change in physics

    #        alpha_tilde : alpha_tilde[q,j] = alpha_hat ** q,q = 0, 1, ... N-1
    #

    # catch cases where N > z.N, need to append extra (N - z.N)
    # copies of any row (all the same) to the remaining rows of alpha_hat
    if z.N < sim_params['N']: # extend the alpha_hat matrix by (N - z.N) rows
        alpha_hat = np.zeros([sim_params['N'], z.CFL.frac[stage, :, :].shape[1]])

        # insert all z.CFL.frac[stage,:,:] rows into alpha_hat;
        alpha_hat[:z.N,:] = z.CFL.frac[stage, :, :]

        # since N > z.N, there are extra rows that need to be filled, store number
        # of extra rows beyond the extent of z.CFL.frac
        N_extra_rows = sim_params['N'] - z.CFL.frac[stage, :, :].shape[0]

        # generate the object of the required size to fit the remaining
        # (N_extra_rows by vz.N) submatrix of alpha_hat that needs to be
        # filled.

        # To do so, compute an outer product with any row of z.CFL.frac[stage,:,:]
        # (all rows i, z.CFL.frac[stage,i,:] are the same as the velocity spans
        # the columns j in z.CFL.frac[stage,i,j] and the CFL numbers depends on
        # the velocity, not on configuration. We copy the first row for
        # obviousness
        alpha_hat_extras = np.outer( np.ones(N_extra_rows), z.CFL.frac[stage, 0, :])
        alpha_hat[z.N:z.N + N_extra_rows, :] = alpha_hat_extras

    else: # z.N >= N, keep N rows for alpha_hat
        alpha_hat = z.CFL.frac[stage, :sim_params['N'], :vz.N]


    # see DECSKS-09 part 1, each vector beta depends on the sign of the CFL:
    #    beta_pos = A_pos.dot(z.CFL.frac[z.CFL.frac >= 0])
    #    beta_neg = A_neg.dot(z.CFL.frac[z.CFL.frac > 0])
    #
    # where A_pos and A_neg are stored arrays constructed in lib.read before
    # the simulation starts:

    # local copies of A matrices
    A_neg = sim_params['A_matrix']['-1']
    A_pos = sim_params['A_matrix']['1']

    # assembling alpha_tilde
    #     alpha_tilde[i,j] = alpha[i,j] ** (q+1) / (q+1)!,
    # where q = 0, 1, ..., N-1
    #
    # we perform the arithmetic operations point-wise via ufuncs
    # by using an array N_arr containg the values of q
    N_arr = np.outer(np.arange(1,sim_params['N']+1), np.ones([1,vz.N]))
    alpha_tilde = ma.array(alpha_hat ** N_arr / scipy.misc.factorial(N_arr))


    # create masks, note we must check the sign on alpha_hat = truncated or
    # expanded z.CFL.frac, we cannot check the sign on alpha_tilde as we have
    # modified the values by raising them to powers, which affects the sign
    mask_neg = (alpha_hat < 0)
    mask_pos = np.logical_not(mask_neg)

    # mask out negative values, leave only positives (>= 0)
    alpha_tilde.mask = mask_neg

    # operate on only positive values (>= 0)
    beta_pos = ma.dot(A_pos, alpha_tilde)

    # mask out positive (>= 0), leave only negatives
    alpha_tilde.mask = mask_pos

    # operate on only negative values
    beta_neg = ma.dot(A_neg, alpha_tilde)

    # consolidate all columns in a single matrix
    B = np.zeros([sim_params['N'], vz.N])

    # wherever beta_neg.mask is True (i.e. marks a negative entry), assign beta_neg, otherwise beta_pos
    B = np.where(mask_neg == True, beta_neg.data, beta_pos.data)

    return B

def correctors(sim_params, z, vz):
    """computes the correction coefficients c for every [i,j]

    inputs:
    sim_params -- (dict) simulation parameters

    outputs:
    c -- (ndarray, ndim=2) correction coefficients
         with shape = (N, z_notadv.N)
         where z_notadv means the not advecting phase
         space variable in the scope of a 2D advection
         implementation

         per DECSKS-09 notation, the tensor c is given by

             c = I_alternating.dot(B)

        where I_alternating.shape = (N,N), and the entries

             I_alternating[i,i] = (-1) ** i
             I_alternating[i,j] = 0 for all i != j

        and the matrix B.shape = (N,vz.N) is the vectors
        of beta correctors (shape = (N,1)) for each value
        of vz.prepointvaluemesh[:,j]
    """
    B = Beta_matrix(sim_params, z, vz)
    c = sim_params['I_alternating'].dot(B)
    return c

def compute_all_correctors_on_a_configuration_variable(sim_params, z, vz, t):
    """computes the correction coefficientxs c for every [i,j]
    in a configuration variable z according to the already
    computed CFL numbers, as stored in the subinstance
    z.CFL

    inputs:
    sim_params -- (dict) simulation parameters
    z -- (instance) configuration variable
    vz -- (instance) physical velocity variable
    t -- (instance) time variable

    outputs:
    c -- (ndarray, ndim=3) each entry c[stage, i, j]
         gives the corrector at each i,j at stage number, stage

         recall from the notebook DECSKS-09, we compute
         the array c per

             c = I_alternating.dot(B)

        where I_alternating.shape = (N,N) and is the identity
        matrix I[q,l] = (-1) ** q if q = l, q = 0, 1, 2, ...
        and is zero otherwise.

        The array B.shape = (N, vz.N) is the assembly of the
        set of N beta correctors (shape = (N,1)) for each value
        of vz.prepointvaluemesh[:,j], i.e. each column of B
        contains the q = 0, 1, ... N-1 corectors, beta, that is
        valid for every i (configuration label) at the given
        j (velocity label)
    """

    # set up dimensions to initialize 3D arrays with
    # recall we have arbitrarily selected split schemes to have
    # the step 'a' to correspond to configuration, and by
    # consequence this means that DECSKS implements 'b' as
    # advecting velocity variables

    dim1 = sim_params['splitting']['number_of_stages']['a'] + 1
    dim2 = sim_params['N']
    dim3 = vz.N

    B = np.zeros((dim2, dim3))
    c = np.zeros((dim1, dim2, dim3))

    for stage in range(1,dim1):
        B = Beta_matrix_for_given_stage(sim_params, z, vz, stage)
        c[stage,:,:] = sim_params['I_alternating'].dot(B)

    return c

def beta(p,a):
    """Returns the pth coefficient beta.

    inputs:
    p -- (int) pth coefficient
    a -- (float) uncorrected fractional distance of MC

    output:

    beta_p -- pth beta coefficient for given vel. or accel.
    """

    beta_p = a ** (p+1) / scipy.misc.factorial(p+1)

    if p > 0:

        for q in range(p):

            if a >= 0:

                beta_p -= beta(q,a) / scipy.misc.factorial(p + 1 - q)

            else: # a < 0

                beta_p -= (-1) ** (p+q) * beta(q,a) / scipy.misc.factorial(p + 1 - q)

    return beta_p
# ........................................................................... #
# DEPRECATED in v2.0

def beta_m(a, B, N):
    A = np.zeros([N,N])
    alpha = np.zeros([N,1])
    for i in range(N):
        alpha[i,0] = a ** (i+1) / scipy.misc.factorial(i+1)
        for j in range(i+1):
            A[i,j] = B[i-j] / scipy.misc.factorial(i-j)
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
