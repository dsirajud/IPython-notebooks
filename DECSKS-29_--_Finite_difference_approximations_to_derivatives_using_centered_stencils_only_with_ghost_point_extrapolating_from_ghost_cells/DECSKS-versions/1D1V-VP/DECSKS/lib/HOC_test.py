def Beta_matrix_v(sim_params, z, vz):
    """constructs the B matrix, whose columns are
    the beta vectors (shape = N x 1) for each value
    of the generalized velocity for the advecting
    variable z.

    See DECSKS-09 part 1 for details of this calculation

    inputs:
    sim_params -- (dict) simulation parameters
    z.CFL.frac -- (ndarray, ndim=2) contains the fractional CFL numbers
                  for every [i,j]
    """

    # This is for v advection only

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

    alpha_hat = z.CFL.frac[:sim_params['N'],:z.CFL.frac.shape[1]]
    alpha_tilde = ma.array(alpha_hat ** N_arr / scipy.misc.factorial(N_arr))

    mask_neg = (alpha_hat < 0)
    mask_pos = np.logical_not(mask_neg)

    # mask out negative values, leave only positives (>= 0)
    alpha_tilde.mask = mask_neg

    # operate on only positive vlaues (>= 0)
    beta_pos = ma.dot(A_pos, alpha_tilde)

    # mask out positive (>= 0), leave only negatives
    alpha_tilde.mask = mask_pos

    # operate on only negative values
    beta_neg = ma.dot(A_neg, alpha_tilde)

    # consolidate all columns in a single matrix
    B = np.zeros([sim_params['N'], sim_params['N' + z.v_str]])

    # wherever beta_neg.mask is False (i.e. unmasked value), assign beta_neg, otherwise beta_pos
    B = np.where(mask_neg == True, beta_neg.data, beta_pos.data)

    return B
