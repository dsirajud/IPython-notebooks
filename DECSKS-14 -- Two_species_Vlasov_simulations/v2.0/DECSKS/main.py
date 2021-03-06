#!/usr/bin/env python
#============================================================================##
#DECSKS - DEterministic Convected Scheme Kinetic Solver for Boltzmann systems #
#-----------------------------------------------------------------------------#
# 1D1V Vlasov-Poisson system, two species                                     #
#                                                                             #
#     __author__     = David Sirajuddin                                       #
#     __version__    = 2.0                                                    #
#     __email__      = sirajuddin@wisc.edu                                    #
#     __status__     = in development                                         #
#                                                                             #
# Python code is crafted with some attention to PEP8 style standards          #
# https://www.python.org/dev/peps/pep-0008/ for                               #
#                                                                             #
#     Python 2.7.3                                                            #
#     NumPy 1.11.0.dev0+fe64f97                                               #
#                                                                             #
#     is not compatible with Python 3+ and/or earlier Numpy releases          #
#                                                                             #
#     coding conventions:                                                     #
#     packages/modules -- lowercase, no underscore                            #
#     classes          -- CapWords                                            #
#     instances        -- lowercase (if no symbol conflict), with underscores #
#     functions        -- lowercase, with underscores unless proper name      #
#                                                                             #
#     other (local) conventions:                                              #
#     numpy arrays     -- lowercase (if no symbol conflict), with underscores #
#     iterators        -- i, j, n used according to f(x_i, v_j, t^n),         #
#     phase space var  -- z = {x, y, z, vx, vy, z}                            #
#=============================================================================#
import _mypath     # adds relative path to sys.path for flexible deployment
import DECSKS
import matplotlib.pyplot as plt
import numpy as np
import time
# =========================================================================== #

#rm_plots = int(raw_input('remove ALL plot files after simulation is done (1 = yes, 0 = no)?: '))
rm_plots = 0
tic = time.clock()

sim_params = DECSKS.lib.read.inputfile('./etc/params_twospecies.dat')

# both species will use same grid x, vx. Can reuse the same vx and ax here
# given serial implementation. In parallel applications, distinct vx_i, vx_e
# ax_i, ax_e may be desirable depending on how the parallelization is approached
x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
vx = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
ax = DECSKS.lib.domain.Setup(sim_params, 'a', 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')

print sim_params['density']
# set up two species
fe, fi = DECSKS.lib.density.setup(sim_params, t, x, vx)

# store total mass for conservation checks, TODO do not need the x phase space variable pass
sim_params['me_0'] = np.sum(np.abs(DECSKS.lib.domain.extract_active_grid(fe[0,:,:], x, sim_params)))
sim_params['mi_0'] = np.sum(np.abs(DECSKS.lib.domain.extract_active_grid(fi[0,:,:], x, sim_params)))

# TODO add this to lib.read, right now you need to make sure you use
# TODO the same mu and tau as in density.setup for sim_params['density']
sim_params['mu'] = 1836.15267389 # mass ratio mi / me for Hydrogen
sim_params['tau'] = 1. # Ti / Te temperature ratio

if sim_params['HOC'] == 'FD':
    sim_params['W'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrices(sim_params,x,vx) # for corrections
    sim_params['W_dn1'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrix_single_derivative(sim_params,x,dn = 1, LTE = 6) # for 6th order Poisson solve

    #DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,0,x,vx)

Plot = DECSKS.lib.plots.PlotSetup(fe, 0, t, x, vx, sim_params, species = 'electron')
Plot(n = 0)
Plot = DECSKS.lib.plots.PlotSetup(fi, 0, t, x, vx, sim_params, species =  'ion')
Plot(n = 0)

    # create and store densities at a chosen x, here we choose x = 0, or x.gridpoints[240 / 2] = x.gridpoints[120]
fe_v = np.zeros([t.N+1, fe.shape[2]])
fi_v = np.zeros([t.N+1, fi.shape[2]])

fe_v[0,:] = fe[0,120,:]
fi_v[0,:] = fi[0,120,:]

print 'simulation has started, status updates are broadcasted after each timestep'

for n in t.stepnumbers:
    print sim_params['split_function_handle']
    fe, fi = eval(sim_params['split_function_handle'])(
        fe, fi,
        t, x, vx, ax,
        n,
        sim_params
        )

    #    Plot = DECSKS.lib.plots.PlotSetup(fe, n, t, x, vx, sim_params, species = 'electron')
    #    Plot(n)
    #    Plot = DECSKS.lib.plots.PlotSetup(fi, n, t, x, vx, sim_params, species =  'ion')
    #    Plot(n)


    # calcs performed and outputs written only if "record outputs? = yes"
    # in ./etc/params.dat
    #    DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,n,x,vx)
    DECSKS.lib.status.check_and_clean(t, n, tic, rm_plots)

    fe_v[n,:] = fe[n,120,:]
    fi_v[n,:] = fi[n,120,:]

for n in range(t.N+1):
    plt.plot(vx.gridvalues,fe_v[n,:], linewidth = 2, label = r'$f_e(t = %g,x=0,v_x)$' % t.times[n])
    plt.grid()
    plt.xlabel(r'$v_x$', fontsize = 18)
    plt.ylabel(r'$f_i(x)$', fontsize = 18)
    plt.legend(loc = 'best')

plt.figure()
for n in range(t.N+1):
    plt.plot(vx.gridvalues,fi_v[n,:], linewidth = 2, label = r'$f_i(t = %g,0,v_x)$' % t.times[n])
    plt.grid()
    plt.xlabel(r'$v_x$', fontsize = 18)
    plt.ylabel(r'$f_i(x)$', fontsize = 18)
    plt.legend(loc = 'best')

plt.show()




toc = time.clock()
simtime = toc - tic
print "simulation completed in %g seconds = %g minutes = %g hours " % (
    simtime,
    simtime/60.,
    simtime/3600.)
# =============================================================================== #
# END
