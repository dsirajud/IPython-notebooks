#!/usr/bin/env python
#============================================================================##
#DECSKS - DEterministic Convected Scheme Kinetic Solver for Boltzmann systems #
#-----------------------------------------------------------------------------#
# CURRENT SETUP: 1D1V, Vlasov-Poisson system with cold ion background         #
#                                                                             #
#     __author__     = David Sirajuddin                                       #
#     __version__    = 1.1                                                    #
#     __email__      = sirajuddin@wisc.edu                                    #
#     __status__     = in development                                         #
#                                                                             #
# Python code is crafted with some attention to PEP8 style standards          #
# for Python 2.7.3                                                            #
#                                                                             #
#     https://www.python.org/dev/peps/pep-0008/                               #
#                                                                             #
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

rm_plots = int(raw_input('remove ALL plot files after simulation is done (1 = yes, 0 = no)?: '))
tic = time.clock()

sim_params = DECSKS.lib.read.inputfile('./etc/params.dat')

x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
v = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')
f = DECSKS.lib.density.setup(sim_params, t, x, v)    # f = f(x_i, v_j, t^n)

# Current case: uniform background (cold) density of ions,
sim_params['ni'] = DECSKS.lib.density.cold_background(f,x,v,sim_params)
# sim_params['W'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrices(sim_params,x,v)
sim_params['m_0'] = np.sum(np.abs((f[0,:x.N,:v.N])))


DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,0,x,v)

#Plot = DECSKS.lib.plots.PlotSetup(f, 0, t, x, v, sim_params)
#Plot(n = 0)

print 'simulation has started, next status update at 10% completion \n'

for n in t.stepnumbers:

    f = DECSKS.lib.split.scheme(
        f,
        t,x,v,
        n,
        sim_params
        )

    #    Plot = DECSKS.lib.plots.PlotSetup(f, n, t, x, v, sim_params)
    #    Plot(n)

    # calcs performed and outputs written only if "record outputs? = yes"
    # in ./etc/params.dat
    #   DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,n,x,v)
    DECSKS.lib.status.check_and_clean(t, n, tic, rm_plots)



toc = time.clock()
print "simulation completed in %g seconds" % (toc - tic)

# =============================================================================== #
# END
