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

f = np.zeros([99,10]) # ad-hoc fix for construction above which does active gridpoints 100 - 1, of no consequence

# n_total = 10, dn = 1, 2, 3, 4, 5, x.N = 99
# f.shape = (x.N, n_total)
# df.shape = (x.N, n_total, dn_total)
# df_approx.shape = W.dot(f).shape = (dn_total, x.N, x.N).dot(x.N, n_total) = (dn_total, x.N, n_total)
# dn = derivative number, n = function label
dn_max = 5
df = np.zeros([dn_max + 1, f.shape[0], f.shape[1]])
df_approx = np.zeros_like(df)
x.prepointvalues = np.linspace(0,2*np.pi, 99)
x.width = 2 * np.pi / 99.

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

import matplotlib.pyplot as plt

W = sim_params['W']
Wx = W['x']

df_approx = Wx.dot(f)

for dn in range(dn_max):
    print "plotting dn - %d, all columns" % dn
    for n in range(6):
        #    plt.plot(x.prepointvalues, f[:,n], linewidth = 2, label = "f, n = %d" % n)
        plt.plot(x.prepointvalues, df_approx[dn, :,n] / x.width ** dn, 'o', linewidth = 2, label = "approx, n = %d, dn = %d" % (n,dn))
        plt.plot(x.prepointvalues, df[dn, :, n], linewidth = 2, label = "approx, n = %d, dn = %d" % (n,dn))
    plt.legend(loc = 'best')
    plt.title('dn = %d' % dn)
    plt.grid()
    plt.figure()

plt.show()




#for n in range(1,f.shape[1]):
#    plt.plot(x.prepointvalues,f[:,n], label = "n = %d" % n)

#plt.legend(loc = 'best')
#plt.title('dn = 1')
#plt.grid()
#plt.figure()
