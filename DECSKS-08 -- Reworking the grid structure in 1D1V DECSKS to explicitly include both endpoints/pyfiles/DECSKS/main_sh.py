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
import sys
import time
# =========================================================================== #
# cwd: /DECSKS
rm_plots = 0 # no input for script shell execution
cleanup_report = 0 # no report of file cleanup printed
tic = time.clock()

params_dir = './etc/NGC_convergence_tests/GB3/F8_2nd/params/'
params_filename = sys.argv[1]
params_filepath = params_dir + params_filename

out_dir = './etc/NGC_convergence_tests/GB3/F8_2nd/outfiles/'

sim_params = DECSKS.lib.read.inputfile(params_filepath)
sim_params['BC'] = 'periodic'
x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')
f = DECSKS.lib.density.setup(sim_params, t, x)    # f = f(t^n, x_i)
v = DECSKS.lib.velocityfields.Velocity('const',x)

if sim_params['HOC'] == 'FD':
    sim_params['W'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrices(sim_params,x)

Plot = DECSKS.lib.plots.PlotSetup(f, 0, t, x, v, sim_params)
Plot(n = 0)

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
print "Nx%d L2 error =" % x.Ngridpoints
print L2
toc = time.clock()
simtime = toc - tic
print "simulation completed in %g seconds = %g minutes = %g hours " % (simtime,
                                                                       simtime/60.,
                                                                       simtime/3600.)
# =============================================================================== #
# END
