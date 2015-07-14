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
import time
# =========================================================================== #

rm_plots = int(raw_input('remove ALL plot files after simulation is done (1 = yes, 0 = no)?: '))
tic = time.clock()

sim_params = DECSKS.lib.read.inputfile('./etc/params.dat')
sim_params['BC'] = 'periodic'
x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
v = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')
f = DECSKS.lib.density.setup(sim_params, t, x, v)    # f = f(x_i, v_j, t^n)

# Current case: uniform background (cold) density of ions,
sim_params['ni'] = DECSKS.lib.density.cold_background(f,x,v,sim_params)

if sim_params['HOC'] == 'FD':
    sim_params['W'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrices(sim_params,x,v)
    sim_params['W_dn1'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrix_single_derivative(sim_params,x,dn = 1, LTE = 6)

DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,0,x,v)

#Plot = DECSKS.lib.plots.PlotSetup(f, 0, t, x, v, sim_params)
#Plot(n = 0)

print 'simulation has started, status updates are broadcasted after each timestep'

for n in t.stepnumbers:

    f = DECSKS.lib.split_poisson.scheme(
        f,
        t,x,v,
        n,
        sim_params
        )

    #    Plot = DECSKS.lib.plots.PlotSetup(f, n, t, x, v, sim_params)
    #    Plot(n)

    # calcs performed and outputs written only if "record outputs? = yes"
    # in ./etc/params.dat
    DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,n,x,v)
    DECSKS.lib.status.check_and_clean(t, n, tic, rm_plots)

toc = time.clock()
simtime = toc - tic
print "simulation completed in %g seconds = %g minutes = %g hours " % (
    simtime,
    simtime/60.,
    simtime/3600.)
# =============================================================================== #
# END
