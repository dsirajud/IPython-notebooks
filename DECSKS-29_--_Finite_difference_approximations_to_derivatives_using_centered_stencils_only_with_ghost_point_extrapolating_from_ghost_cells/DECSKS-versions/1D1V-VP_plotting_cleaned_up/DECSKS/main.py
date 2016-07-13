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
#                      -- k is used for contiguous postpointsm k = 0,1
#     phase space var  -- z = {x, y, z, vx, vy, z}                            #
#=============================================================================#
import _mypath     # adds relative path to sys.path for flexible deployment
import DECSKS
import extra_plots
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys
import time
# =========================================================================== #

#rm_plots = int(raw_input('remove ALL plot files after simulation is done (1 = yes, 0 = no)?: '))
rm_plots = 0
tic = time.clock()


sim_name = 's17-04'
sim_params = DECSKS.lib.read.inputfile('./etc/params_' + sim_name + '.dat')

# both species will use same grid x, vx. Can reuse the same vx and ax here
# given serial implementation. In parallel applications, distinct vx_i, vx_e
# ax_i, ax_e may be desirable depending on how the parallelization is approached
x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
vx = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
ax = DECSKS.lib.domain.Setup(sim_params, 'a', 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')

# set up two species
fe, fi = DECSKS.lib.density.setup(sim_params, t, x, vx) # NOTE mu and tau in ion density must match those just below

print sim_params['compute_electric_potential_phi_handle'][x.str]
print sim_params['phi_BC']['x']
print "CHECK BOUNDARIES MATCH INPUTS"
print "BC type on f on x"
print sim_params['BC']['f']['x']['type']
print "BC type on f on vx"
print sim_params['BC']['f']['vx']['type']
print "left BC on f on x"
print sim_params['BC']['f']['x']['lower']
print "right BC on f on x"
print sim_params['BC']['f']['x']['upper']
print "left BC on f on vx"
print sim_params['BC']['f']['vx']['lower']
print "right BC son f on vx"
print sim_params['BC']['f']['vx']['upper']

if sim_params['BC']['f']['x']['upper'] == 'source':
    DECSKS.lib.sources.setup(x, vx, t, sim_params)

    #print "plotting the in-flow sources"
    #Plot = DECSKS.lib.plots.PlotSetup(sim_params['sources']['f_Se'][1,:,:], 0, t, x, vx, sim_name, sim_params, species = 'source electrons from right')
    #Plot(n = 0)

# Create plotting class objects (contains __call__ plot routine and __init__ sets up the file name and save location)
PlotElectron = DECSKS.lib.plots.PlotSetup(fe, t, x, vx, sim_name, sim_params, species = 'electron', quantity = 'distribution_function')
PlotIon = DECSKS.lib.plots.PlotSetup(fi, t, x, vx, sim_name, sim_params, species =  'ion', quantity = 'distribution_function')

# plot and save
PlotElectron(fe, n = 0, sim_name = sim_name)
PlotIon(fi, n = 0, sim_name = sim_name)

phi_2D = eval(sim_params['compute_electric_potential_phi_handle'][x.str])(fe, fi, x, vx, 0, sim_params)
phi = phi_2D[:,0] # all columns are identical

# create plot object and call to plot/save
PlotPhi = DECSKS.lib.plots.PlotSetup(phi, t, x, sim_name = sim_name, sim_params = sim_params, quantity = 'potential')
PlotPhi(phi, n = 0, sim_name = sim_name)

#DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,fe, fi, 0, x, vx, sim_params['mu'])

IQ = np.sum(-fe + fi) * vx.width * x.width
print "\n initial charge in system"
print IQ

print 'simulation has started, status updates are broadcasted after each timestep'
for n in t.stepnumbers:
    fe, fi = DECSKS.lib.split.scheme(
        fe, fi,
        t, x, vx, ax,
        n,
        sim_params
        )

    PlotElectron(fe, n, sim_name = sim_name) # plot fe at time step n
    PlotIon(fi, n, sim_name = sim_name)      # plot fi at time step n

    # calc and plot/save phi
    phi_2D = eval(sim_params['compute_electric_potential_phi_handle'][x.str])(fe, fi, x, vx, n, sim_params)
    phi = phi_2D[:,0] # all columns are identical
    PlotPhi(phi, n, sim_name)

    #-----------------------------------------------------------------#
    # plot densities a few times throughout the simulation
    #-----------------------------------------------------------------#
    if n == t.N/1000:

        extra_plots.densities(fe,fi, x, vx, t, n, sim_name, sim_params)

    elif n == t.N/10:

        extra_plots.densities(fe,fi, x, vx, t, n ,sim_name, sim_params)

    elif n == t.N/4:

        extra_plots.densities(fe,fi, x, vx, t, n , sim_name,sim_params)

    elif n == t.N/2:

        extra_plots.densities(fe,fi, x, vx, t, n , sim_name, sim_params)

    elif n == 3*t.N/4:

        extra_plots.densities(fe,fi, x, vx, t, n , sim_name,sim_params)

    elif n == 9*t.N/10:

        extra_plots.densities(fe,fi, x, vx, t, n , sim_name,sim_params)

    elif n == .95*t.N:

        extra_plots.densities(fe,fi, x, vx, t, n , sim_name,sim_params)

    elif n == t.N:

        extra_plots.densities(fe,fi, x, vx, t, n , sim_name,sim_params)
    #-----------------------------------------------------------------#

    # calcs performed and outputs written only if "record outputs? = yes"
    # in ./etc/params.dat
    #    DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,fe, fi, n, x, vx, sim_params['mu'])
    DECSKS.lib.status.check_and_clean(t, n, tic, rm_plots)

toc = time.clock()
simtime = toc - tic
print "simulation completed in %g seconds = %g minutes = %g hours " % (
    simtime,
    simtime/60.,
    simtime/3600.)
# =============================================================================== #
# END
