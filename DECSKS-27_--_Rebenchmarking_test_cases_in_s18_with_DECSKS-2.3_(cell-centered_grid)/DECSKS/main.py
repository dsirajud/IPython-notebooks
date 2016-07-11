#!/usr/bin/env python
#============================================================================##
#DECSKS - DEterministic Convected Scheme Kinetic Solver for Boltzmann systems #
#-----------------------------------------------------------------------------#
# 1D1V Vlasov-Poisson system, two species                                     #
#                                                                             #
#     __author__     = David Sirajuddin                                       #
#     __version__    = 2.3                                                    #
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
import matplotlib
import numpy as np
import time
# =========================================================================== #

#rm_plots = int(raw_input('remove ALL plot files after simulation is done (1 = yes, 0 = no)?: '))
rm_plots = 0
tic = time.clock()

sim_params = DECSKS.lib.read.inputfile('./etc/params_s07-01.dat')

# both species will use same grid x, vx. Can reuse the same vx and ax here
# given serial implementation. In parallel applications, distinct vx_i, vx_e
# ax_i, ax_e may be desirable depending on how the parallelization is approached
x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
vx = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
ax = DECSKS.lib.domain.Setup(sim_params, 'a', 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')

print x.prepointvalues
print x.width
print x.prepointvalues[0] - x.width / 2
print x.prepointvalues[-1] + x.width / 2
print x.N

print vx.prepointvalues
print vx.width
print vx.prepointvalues[0] - vx.width / 2
print vx.prepointvalues[-1] + vx.width / 2
print vx.N


# set up two species
fe, fi = DECSKS.lib.density.setup(sim_params, t, x, vx) # NOTE mu and tau in ion density must match those just below


ne_avg = np.sum(fe)*x.width * vx.width / x.L
print ne_avg

# store total mass for conservation checks, TODO do not need the x phase space variable pass in this function
sim_params['me_0'] = np.sum(fe)
sim_params['mi_0'] = np.sum(fi)

print "TIME ZERO, masses are"
print "fe = %g" % np.sum(fe)
print "fi = %g" % np.sum(fi)

print x.prepointvaluemesh
print x.prepointmesh
print x.width == x.L / x.N
print vx.prepointvaluemesh
print vx.prepointmesh
print vx.width == vx.L / vx.N
# CAUTION: make sure you use the same mu and tau in lib.density for fe, fi as indicated in params.dat input deck
# as the mu specified in params.dat will be used to compute the ion acceleration ax ~ 1 / mu * Ex compared to the
# electron acceleration term ax ~ -Ex

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

#Ex_2D = eval(sim_params['compute_electric_field_orchestrator_handle']['x'])(fe, fi, x, vx, 0, sim_params)
#Ex = Ex_2D[:,0] # all columns of the original 2D array are the same, choose the zeroeth column arbitrarily

#phi = eval(sim_params['compute_electric_potential_phi_handle'][x.str])(fe, fi, x, vx, 0, sim_params)
#phi = phi[:,0]
#print phi.shape
DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,fe, fi, 0, x, vx, sim_params['mu'])

#matplotlib.pyplot.plot(x.gridvalues, phi, linewidth = 2, color = 'blue')
#matplotlib.pyplot.grid()
##matplotlib.pyplot.axis([x.gridvalues[0], x.gridvalues[-1], 0, 35])
#matplotlib.pyplot.xlabel(r'position $x$', fontsize = 18)
#matplotlib.pyplot.ylabel(r'$\phi (t^n,x)$', fontsize = 18)
#matplotlib.pyplot.title(r's18-21d Electric potential $\phi (x)$: $N_x$ = %d, $N_v$ = %d, $t^n$ = %2.3f, n = %03d' % (sim_params['active_dims'][0], sim_params['active_dims'][1], 0.0, 0))
#it_str = 'it%05d' % 0
#matplotlib.pyplot.savefig('./plots/' + 'phi_s18-21_' + it_str)
#matplotlib.pyplot.clf()

#matplotlib.pyplot.plot(x.gridvalues, Ex, linewidth = 2, color = 'blue')
#matplotlib.pyplot.grid()
#matplotlib.pyplot.axis([x.gridvalues[0], x.gridvalues[-1], 0, 35])
#matplotlib.pyplot.xlabel(r'position $x$', fontsize = 18)
#matplotlib.pyplot.ylabel(r'$E (t^n,x)$', fontsize = 18)
#matplotlib.pyplot.title(r's18-21d Electric field $E (x)$: $N_x$ = %d, $N_v$ = %d, $t^n$ = %2.3f, n = %03d' % (sim_params['active_dims'][0], sim_params['active_dims'][1], 0.0, 0))
#it_str = 'it%05d' % 0
#matplotlib.pyplot.savefig('./plots/' + 'Ex_s18-21_' + it_str)
#matplotlib.pyplot.clf()

#print sim_params['sigma']['x']['lower']
#print sim_params['sigma']['x']['upper']

print 'simulation has started, status updates are broadcasted after each timestep'

for n in t.stepnumbers:
    fe, fi = DECSKS.lib.split.scheme(
        fe, fi,
        t, x, vx, ax,
        n,
        sim_params
        )

    #    sim_params['sigma_n']['x']['lower'][n] = sim_params['sigma']['x']['lower']
    #    sim_params['sigma_n']['x']['upper'][n] = sim_params['sigma']['x']['upper']

    #    print sim_params['sigma']['x']['lower']
    #    print sim_params['sigma']['x']['upper']

    #    Plot = DECSKS.lib.plots.PlotSetup(fe, n, t, x, vx, sim_params, species = 'electron')
    #    Plot(n)
    #    Plot = DECSKS.lib.plots.PlotSetup(fi, n, t, x, vx, sim_params, species =  'ion')
    #    Plot(n)

    #    Ex_2D = eval(sim_params['compute_electric_field_orchestrator_handle']['x'])(fe, fi, x, vx, n, sim_params)
    #    Ex = Ex_2D[:,0] # all columns of the original 2D array are the same, choose the zeroeth column arbitrarily

    #    phi = eval(sim_params['compute_electric_potential_phi_handle'][x.str])(fe, fi, x, vx, 0, sim_params)
    #    phi = phi[:,0]

    #    matplotlib.pyplot.plot(x.gridvalues, phi, linewidth = 2, color = 'blue')
    #    matplotlib.pyplot.grid()
    #    #    matplotlib.pyplot.axis([x.gridvalues[0], x.gridvalues[-1], 0, 35])
    #    matplotlib.pyplot.xlabel(r'position $x$', fontsize = 18)
    #    matplotlib.pyplot.ylabel(r'$\phi (t^n,x)$', fontsize = 18)
    #    matplotlib.pyplot.title(r's18-21d Electric potential $\phi (t^n, x)$: $N_x$ = %d, $N_v$ = %d, $t^n$ = %2.3f, n = %03d' % (sim_params['active_dims'][0], sim_params['active_dims'][1], n*t.width, n))
    #    it_str = 'it%05d' % n
    #    matplotlib.pyplot.savefig('./plots/' + 'phi_s18-21_' + it_str)
    #    matplotlib.pyplot.clf()


    #    matplotlib.pyplot.plot(x.gridvalues, Ex, linewidth = 2, color = 'blue')
    #    matplotlib.pyplot.grid()
    #    #    matplotlib.pyplot.axis([x.gridvalues[0], x.gridvalues[-1], 0, 35])
    #    matplotlib.pyplot.xlabel(r'position $x$', fontsize = 18)
    #    matplotlib.pyplot.ylabel(r'$E (t^n,x)$', fontsize = 18)
    #    matplotlib.pyplot.title(r's18-21d Electric field $E (x)$: $N_x$ = %d, $N_v$ = %d, $t^n$ = %2.3f, n = %03d' % (sim_params['active_dims'][0], sim_params['active_dims'][1], n*t.width, n))
    #    it_str = 'it%05d' % n
    #    matplotlib.pyplot.savefig('./plots/' + 'Ex_s18-21_' + it_str)
    #    matplotlib.pyplot.clf()

    # calcs performed and outputs written only if "record outputs? = yes"
    # in ./etc/params.dat
    DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,fe, fi, n, x, vx, sim_params['mu'])
    DECSKS.lib.status.check_and_clean(t, n, tic, rm_plots)


    #sigma_n_left = sim_params['sigma_n']['x']['lower']
    #sigma_n_right = sim_params['sigma_n']['x']['upper']

    #plt.plot(t.times, sigma_n_left, linewidth = 2, label = r'$\sigma (t, x= -10)$')
    #plt.plot(t.times,sigma_n_right, linewidth = 2, label = r'$\sigma (t, x= +10)$')
    #plt.grid()
    #plt.xlabel(r'time step $n$', fontsize = 18)
    #plt.ylabel(r'$\sigma (t,x)$', fontsize = 18)
    #plt.legend(loc = 'best')

#phi_left = sim_params['sigma_n']['x']['lower'] # E = -1/2 sigma, phi = 1/2 sigma, here sigma = ni - ne
#phi_right = sim_params['sigma_n']['x']['upper']
#plt.plot(trange,phi_left, linewidth = 2, label = r'$\phi (t, x= -10)$')
#plt.plot(trange,phi_right, linewidth = 2, label = r'$\phi (t, x= +10)$')
#plt.grid()
#plt.xlabel(r'time step $n$', fontsize = 18)
#plt.ylabel(r'$\phi (t,x)$', fontsize = 18)
#plt.legend(loc = 'best')

#plt.show()

toc = time.clock()
simtime = toc - tic
print "simulation completed in %g seconds = %g minutes = %g hours " % (
    simtime,
    simtime/60.,
    simtime/3600.)
# =============================================================================== #
# END
