import numpy as np
import numpy.ma as ma
import scipy.misc

N = 21
Nx = 24
Nvx = 49


N_arr = np.outer(np.arange(1,N+1), np.ones([1,Nx]))

alpha = np.random.randn(Nx, Nvx)
alpha = np.transpose(alpha) # v advection

alpha_hat = alpha[:N, : alpha.shape[1]]

alpha_tilde = ma.array(alpha_hat ** N_arr / scipy.misc.factorial(N_arr))
