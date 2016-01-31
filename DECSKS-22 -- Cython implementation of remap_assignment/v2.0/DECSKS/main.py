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

sim_params = DECSKS.lib.read.inputfile('./etc/params.dat')

x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
vx = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
ax = DECSKS.lib.domain.Setup(sim_params, 'a', 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')

# set up two species
fe, fi = DECSKS.lib.density.setup(sim_params, t, x, vx) # NOTE mu and tau in ion density must match those just below

# store total mass for conservation checks, TODO do not need the x phase space variable pass in this function
sim_params['me_0'] = np.sum(np.abs(DECSKS.lib.domain.extract_active_grid(fe, sim_params)))
sim_params['mi_0'] = np.sum(np.abs(DECSKS.lib.domain.extract_active_grid(fi, sim_params)))

# TODO add this to lib.read, right now you need to make sure you use
# TODO the same mu and tau as in density.setup for sim_params['density']
sim_params['mu'] = 1836.15267389 # mass ratio mi / me for Hydrogen, needed for ion acceleration term in lib.split

    #DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,0,x,vx)

Plot = DECSKS.lib.plots.PlotSetup(fe, 0, t, x, vx, sim_params, species = 'electron')
Plot(n = 0)
Plot = DECSKS.lib.plots.PlotSetup(fi, 0, t, x, vx, sim_params, species =  'ion')
Plot(n = 0)

print 'simulation has started, status updates are broadcasted after each timestep'
for n in t.stepnumbers:
    fe, fi = DECSKS.lib.split.scheme(
        fe, fi,
        t, x, vx, ax,
        n,
        sim_params
        )


    Plot = DECSKS.lib.plots.PlotSetup(fe, n, t, x, vx, sim_params, species = 'electron')
    Plot(n)
    Plot = DECSKS.lib.plots.PlotSetup(fi, n, t, x, vx, sim_params, species =  'ion')
    Plot(n)

    # calcs performed and outputs written only if "record outputs? = yes"
    # in ./etc/params.dat
    #    DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,n,x,vx)
    DECSKS.lib.status.check_and_clean(t, n, tic, rm_plots)
toc = time.clock()
simtime = toc - tic
print "simulation completed in %g seconds = %g minutes = %g hours " % (
    simtime,
    simtime/60.,
    simtime/3600.)
# =============================================================================== #
# END
