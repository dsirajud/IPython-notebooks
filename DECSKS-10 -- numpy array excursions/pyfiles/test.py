import timeit
import numpy as np
import scipy
from math import factorial

Nx = 768
Nv = 1536
N = 21

# create an example alpha, note alpha[:,j] = const for each j
alpha = np.zeros([Nx,Nv])
for j in range(Nv):
    alpha[:,j] = np.random.rand()*np.ones(Nx)

# create a A_pos and A_neg matrice, which require Bernoulli numbers Bern[p], need p = 0, 1, ..., N - 1

Bern = np.zeros(21)
Bern[0] = 1.0
Bern[1] = -1/2.
Bern[2] = 1/6.
Bern[3] = 0
Bern[4] = -1/30.
Bern[5] = 0
Bern[6] = 1/42.
Bern[7] = 0
Bern[8] = -1/30.
Bern[9] = 0
Bern[10] = 5/66.
Bern[11] = 0
Bern[12] = -691/2730.
Bern[13] = 0
Bern[14] = 7/6.
Bern[15] = 0
Bern[16] = -3617/510.
Bern[17] = 0
Bern[18] = 43867/798.
Bern[19] = 0
Bern[20] = -174611/330.

# Assemble A matrix

A_pos, A_neg = np.zeros([N,N]), np.zeros([N,N])

for i in range(N):
    for j in range(i+1):
        A_pos[i,j] = Bern[i-j]/factorial(i-j) # Bern[k] is the kth Bernoulli number
        if (i - j) == 1:
            A_neg[i,j] = -A_pos[i,j]
        else:
            A_neg[i,j] = A_pos[i,j]

# Create a selector based on the sign of alpha, which is the key to a subdictionary of the sim_params dict

sim_params = {}
sim_params['A_matrix'] = {}
sim_params['A_matrix']['1'] = A_pos
sim_params['A_matrix']['0'] = A_pos
sim_params['A_matrix']['-1'] = A_neg

Bx = np.zeros([N, Nv])
alpha_tilde = np.zeros([N, Nv])

for timestep in range(360):
    for j in range(Nv):
        alpha_tilde[:,j] = alpha[0,j] ** np.arange(1, N + 1) \
        / scipy.misc.factorial(np.arange(1,N + 1))

        # A_matrix.shape = (N,N), alpha_tilde[:,j].shape = (N,1), Bx[:,j].shape = (N,1)
        # and Bx[:,:].shape = (N, v.N)
        Bx[:,j] = sim_params['A_matrix'][str(int(np.sign(alpha[0,j])))].dot(alpha_tilde[:,j])
