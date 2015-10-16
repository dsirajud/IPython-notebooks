import numpy as np

# run this in ipython from DECSKS/
import _mypath     # adds relative path to sys.path for flexible deployment
import DECSKS
import numpy as np
# =========================================================================== #

#rm_plots = int(raw_input('remove ALL plot files after simulation is done (1 = yes, 0 = no)?: '))
#tic = time.clock()

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

    Nx = 768 - 1 # this is an ad-hoc fix to the fact that we are reading params.dat and it has PBCs, hence Nx -= 1 in the end
f = np.zeros([Nx,10]) 
df = np.zeros_like(f)
x.prepointvalues = np.linspace(-2*np.pi / .3,2*np.pi / .3, Nx)
x.width = 1 / .3 * 4 * np.pi / Nx

for n in range(f.shape[1]):
    f[:,n] = np.sin(x.prepointvalues) + n*x.prepointvalues
    df[:,n] = np.cos(x.prepointvalues) + n

import matplotlib.pyplot as plt

W = sim_params['W']
Wx = W['x']
df_approx = Wx[1,:,:].dot(f)
df_approx /= x.width 

for n in range(f.shape[1]):
    plt.plot(x.prepointvalues,df_approx[:,n], 'o', linewidth = 2, label = "approx, n = %d" % n)
    plt.plot(x.prepointvalues, df[:,n], linewidth = 2, label = 'n = %d' % n)

plt.legend(loc = 'best')
plt.title('dn = 1')
plt.grid()
