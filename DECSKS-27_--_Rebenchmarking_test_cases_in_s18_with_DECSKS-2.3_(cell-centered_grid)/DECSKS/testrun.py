import _mypath     # adds relative path to sys.path for flexible deployment
import DECSKS
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import time
sim_params = DECSKS.lib.read.inputfile('./etc/params_s18-20.dat')

# both species will use same grid x, vx. Can reuse the same vx and ax here
# given serial implementation. In parallel applications, distinct vx_i, vx_e
# ax_i, ax_e may be desirable depending on how the parallelization is approached
x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
vx = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
ax = DECSKS.lib.domain.Setup(sim_params, 'a', 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')

fe, fi = DECSKS.lib.density.setup(sim_params, t, x, vx) # NOTE mu and tau in ion density must match those just below

ne_avg = np.sum(fe)*x.width * vx.width / x.L
print ne_avg

# store total mass for conservation checks, TODO do not need the x phase space variable pass in this function
sim_params['me_0'] = np.sum(fe)
sim_params['mi_0'] = np.sum(fi)

print "TIME ZERO, masses are"
print "fe = %g" % np.sum(fe)
print "fi = %g" % np.sum(fi)
print sim_params['BC']['f']['x']['type']
print sim_params['BC']['f']['vx']['type']
print sim_params['BC']['f']['x']['lower']
print sim_params['BC']['f']['x']['upper']
print sim_params['BC']['f']['vx']['lower']
print sim_params['BC']['f']['vx']['upper']
#print sim_params['compute_electric_potential_phi_handle'][x.str] # = None if fourier solver
# print sim_params['phi_BC']['x'] # = None if fourier solver

Plot = DECSKS.lib.plots.PlotSetup(fe, 0, t, x, vx, sim_params, species = 'electron')
Plot(n = 0)
Plot = DECSKS.lib.plots.PlotSetup(fi, 0, t, x, vx, sim_params, species =  'ion')
Plot(n = 0)

print "integral (fe - fi) dx dvx = %g" % (np.sum(-fe + fi) * x.width * vx.width)

print "integral 2 (f dx dvx - fi dx dvx) = %g" % ((np.sum(-fe) * x.width * vx.width + np.sum(fi) * x.width * vx.width))





