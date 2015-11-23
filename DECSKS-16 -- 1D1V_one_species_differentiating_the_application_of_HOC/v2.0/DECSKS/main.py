#!/usr/bin/env python
#============================================================================##
#DECSKS - DEterministic Convected Scheme Kinetic Solver for Boltzmann systems #
#-----------------------------------------------------------------------------#
# CURRENT SETUP: 1D1V, Vlasov-Poisson system with cold ion background         #
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
f = DECSKS.lib.density.setup(sim_params, t, x, vx)    # f = f(x_i, v_j, t^n)

sim_params['m_0'] = np.sum(np.abs(DECSKS.lib.domain.extract_active_grid(f[0,:,:], x, sim_params))) # initial total mass
sim_params['ni'] = DECSKS.lib.density.cold_background(f,x,vx,sim_params) # cold background case

DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,0,x,vx)
Plot = DECSKS.lib.plots.PlotSetup(f, 0, t, x, vx, sim_params) # title, axes labels, etc.
Plot(n = 0)

print '\nstatus updates are broadcasted after each time step\n'

for n in t.stepnumbers:
    f = DECSKS.lib.split.scheme(
        f,
        t,
        x, vx, ax,
        n,
        sim_params
        )

    Plot = DECSKS.lib.plots.PlotSetup(f, n, t, x, vx, sim_params) # title, axes labels, etc.
    Plot(n)

    # calcs performed and outputs written only if "record outputs? = yes"
    # in ./etc/params.dat
    DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,n,x,vx)
    DECSKS.lib.status.check_and_clean(t, n, tic, rm_plots)

toc = time.clock()
simtime = toc - tic
print "simulation completed in %g seconds = %g minutes = %g hours " % (
    simtime,
    simtime/60.,
    simtime/3600.)
# =============================================================================== #
# END
