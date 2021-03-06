# run this in ipython from DECSKS/

import numpy as np
import _mypath     # adds relative path to sys.path for flexible deployment
import DECSKS
import time

sim_params = DECSKS.lib.read.inputfile('./etc/params.dat')

x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
vx = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')

# store total mass for conservation checks

z = x
vz = vx

# The following W and W_dn1 excerpt is from main.py
if sim_params['HOC'] == 'FD':
    sim_params['W'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrices(sim_params,x,vx)
    sim_params['W_dn1'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrix_single_derivative(sim_params,x,dn = 1, LTE = 6)


Nx = 768 - 1
f = np.zeros([Nx,10]) # ad-hoc fix for construction above which does active gridpoints 100 - 1, of no consequence

# n_total = 10, dn = 1, 2, 3, 4, 5, x.N = 99
# f.shape = (x.N, n_total)
# df.shape = (x.N, n_total, dn_total)
# df_approx.shape = W.dot(f).shape = (dn_total, x.N, x.N).dot(x.N, n_total) = (dn_total, x.N, n_total)
# dn = derivative number, n = function label
dn_max = 5
df = np.zeros([dn_max + 1, f.shape[0], f.shape[1]])
x.prepointvalues = np.linspace(-2 * np.pi / .3, 2*np.pi / .3, Nx)
x.width = 1 / .3 * 4 * np.pi / Nx

dn = 0
for n in range(f.shape[1]):
    f[:,n] = np.sin(n * x.prepointvalues)

dn = 1
for n in range(f.shape[1]):
    df[dn, :, n] = n * np.cos(n * x.prepointvalues)

dn = 2
for n in range(f.shape[1]):
    df[dn, :,n] = -n **2 * np.sin(n * x.prepointvalues)

dn = 3
for n in range(f.shape[1]):
    df[dn, :,n] = - n ** 3 *  np.cos(n * x.prepointvalues)

dn = 4
for n in range(f.shape[1]):
    df[dn, :,n] = n ** 4 * np.sin(n * x.prepointvalues)

dn = 5
for n in range(f.shape[1]):
    df[dn, :,n] = n ** 5 * np.cos(n * x.prepointvalues)

W = sim_params['W']
Wx = W['x']


t = time.time()
df_approx = np.zeros([dn_max, f.shape[0], f.shape[1]])
for dn in range(1,dn_max):
    df_approx[dn,:,:] = Wx[dn,:,:].dot(f)
td = time.time() - t

print "looping through each dn + dot product 2D matrices completed in %.3f s" % td

t = time.time()
df_approx2 = Wx.dot(f)
td = time.time() - t

print "a single 3D dot product completed in %.3f s" % td



# timeit calls

# looping version

#import timeit
#stmt = "derivative_loop(Wx, f, 5)"
#setup = "from derivative_timeit import setup, derivative_loop, derivative_single_dot; Wx, f = setup()"
#print(timeit.timeit(stmt, setup, number = 1))

# single dot version

#import timeit
#stmt = "derivative_single_dot(Wx, f, 5)"
#setup = "from derivative_timeit import setup, derivative_loop, derivative_single_dot; Wx, f = setup()"
#print(timeit.timeit(stmt, setup, number = 1))
