#!/usr/bin/env python
#============================================================================##
#DECSKS - DEterministic Convected Scheme Kinetic Solver for Boltzmann systems #
#-----------------------------------------------------------------------------#
# CURRENT SETUP: 1D1V, Vlasov-Poisson system with cold ion background         #
#                                                                             #
#     __author__     = David Sirajuddin                                       #
#     __version__    = 1.2                                                    #
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
cleanup_report = int(raw_input('provide full report of all removed files (1 = yes, 0 = no)?: '))
tic = time.clock()

sim_params = DECSKS.lib.read.inputfile('./etc/params.dat')
sim_params['BC'] = 'periodic'
x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')
f = DECSKS.lib.density.setup(sim_params, t, x)    # f = f(t^n, x_i)
v = DECSKS.lib.velocityfields.Velocity('const',x)

print x.gridvalues
print "this is main.py, not main_sh.py"
print "x.N = %d" % x.N
print "x.Ngridpoints = %d" % x.Ngridpoints

sim_params['W'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrices(sim_params,x)

#DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,0,x,t)

Plot = DECSKS.lib.plots.PlotSetup(f, 0, t, x, v, sim_params)
Plot(n = 0)
print Plot
print 'simulation has started, next status update at 10% completion \n'

for n in t.stepnumbers:

    x.MCs   = x.generate_Lagrangian_mesh(x.prepointvalues, v.prepointvalues, t.width)
    f[n,:] = DECSKS.lib.convect.scheme(
        f[n-1,:],
        x, n,
        sim_params
        )

    Plot = DECSKS.lib.plots.PlotSetup(f, n, t, x, v, sim_params)
    Plot(n)

    # calcs performed and outputs written only if "record outputs? = yes"
    # in ./etc/params.dat
    #    DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,n,x)
    DECSKS.lib.status.check_and_clean(t, n, tic, rm_plots, cleanup_report)

L2 = DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,t.N,x,t)
print "L2 error ="
print L2
toc = time.clock()
print "simulation completed in %g seconds" % (toc - tic)

# =============================================================================== #
# END
