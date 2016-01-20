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

    # operate on only positive vlaues (>= 0)
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

def Beta_matrix_3D(sim_params, z, vz, s):
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
    if z.CFL.frac[s,:,:].shape[0] < sim_params['N']:
        alpha_hat = np.zeros([sim_params['N'], z.CFL.frac[s,:,:].shape[1]])
        alpha_hat[:z.CFL.frac[s,:,:].shape[0],:] = z.CFL.frac[s,:,:]
        N_extra_rows = sim_params['N'] - z.CFL.frac[s,:,:].shape[0]

        # generate the object of the required size to fit the remaining submatrix
        # in alpha_hat not already filled, compute an outer product with any row (all same)
        # here, we take the first row for obviousness
        alpha_hat_extras = np.outer( np.ones(N_extra_rows), z.CFL.frac[s,0,:])
        alpha_hat[z.CFL.frac[s,:,:].shape[0]:z.CFL.frac[s,:,:].shape[0] + N_extra_rows,:] = alpha_hat_extras

    else:
        alpha_hat = z.CFL.frac[s,:sim_params['N'],:z.CFL.frac[s,:,:].shape[1]]

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
    B = np.zeros([sim_params['N'], vz.N])

    # wherever beta_neg.mask is False (i.e. unmasked value), assign beta_neg, otherwise beta_pos
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

def correctors_on_configuration(sim_params, z, vz, t):
    """computes the correction coefficients c for every [i,j]
    in a configuration variable z according to velocity vz

    inputs:
    sim_params -- (dict) simulation parameters

    outputs:
    c -- (ndarray, ndim=2) correction coefficients
         with shape = (N, vz.N)

         per DECSKS-09 notation, the tensor c is given by

             c = I_alternating.dot(B)

        where I_alternating.shape = (N,N), and the entries

             I_alternating[i,i] = (-1) ** i
             I_alternating[i,j] = 0 for all i != j

        and the matrix B.shape = (N,vz.N) is the vectors
        of beta correctors (shape = (N,1)) for each value
        of vz.prepointvaluemesh[:,j]
    """

    # compute CFL numbers for every substage in a split scheme

    dim1 = sim_params['splitting']['number_of_substeps']['a'] + 1 # index = 0 not used
    dim2 = sim_params['N']
    dim3 = vz.N

    B = np.zeros([dim1, dim2, dim3])
    c = np.zeros([dim1, dim2, dim3])
    
    for s in range(1,dim1):
        B[s,:,:] = Beta_matrix_3D(sim_params, z, vz, s)
        c[s,:,:] = sim_params['I_alternating'].dot(B[s,:,:])

    return c
