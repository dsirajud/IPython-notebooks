#!/usr/bin/env python
#=============================================================================#
#DECSKS - DEterministic Convected Scheme Kinetic Solver for Boltzmann systems #
#-----------------------------------------------------------------------------#
# 1D1V Vlasov-Poisson system, two species                                     #
#                                                                             #
#     __author__     = David Sirajuddin                                       #
#     __version__    = 2.4                                                    #
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
#                      -- k is used for contiguous postpointsm k = 0,1        #
#     phase space var  -- z = {x, y, z, vx, vy, z}                            #
#=============================================================================#
import _mypath     # adds relative path to sys.path for flexible deployment
import DECSKS
import notifications
import numpy as np
import time
# =========================================================================== #
tic = time.clock() # start global clock (tracks full simulation processor time)
# *************************************************************************** #
# (0) SIMULATION SPECIFICATION                                                #
# --------------------------------------------------------------------------- #
# Specify simulation and grab inputs indicated in etc/params_sim_name.dat
sim_name = 's17-04'
sim_params = DECSKS.lib.read.inputfile(sim_name)
# --------------------------------------------------------------------------- #
# *************************************************************************** #
# (1) CREATE GRID AND INITIALIZE DISTRIBUTIONS                                #
# --------------------------------------------------------------------------- #
# instantiate phase space variables (x, vx), acceleration (ax), and time (t)
x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
vx = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
ax = DECSKS.lib.domain.Setup(sim_params, 'a', 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')

# initialize distribution functions
fe, fi = DECSKS.lib.density.setup(sim_params, t, x, vx) # NOTE mu and tau in ion density must match those just below

# print out print some parameters of interest if unit testing
notifications.announce(fe, fi, x, vx, sim_params, print_to_screen = True)
# --------------------------------------------------------------------------- #
# *************************************************************************** #
# (2) MISCELLANEOUS AD-HOC ACTIONS (IF ANY)
# --------------------------------------------------------------------------- #
# This section is used to rapidly introduce new features for unit testing
# and/or results, i.e. functionality that has not been generalized that
# may pertain only to as few as one simulation

# add source setup to lib.read later
if sim_params['BC']['f']['x']['upper'] == 'source':
    DECSKS.lib.sources.setup(x, vx, t, sim_params)
    #print "plotting the in-flow sources"

#fe_ext, fi_ext = ghostgrid.setup(fe, fi, x, vx)
# --------------------------------------------------------------------------- #
# *************************************************************************** #
# (3) WRITE QUANTITY BIOGRAPHIES AND DATA
# --------------------------------------------------------------------------- #
# write fe, fi to files (./etc/outputs/sim_name/)
DECSKS.lib.write.tofiles(sim_name, n = 0,
                         fe = fe,
                         fi = fi
                         )

# write grid to files
DECSKS.lib.write.tofiles(sim_name,
                         xgridvalues = x.gridvalues,
                         vxgridvalues = vx.gridvalues,
                         times = t.times
                         )

if sim_params['HOC']['x'] == 'FD':
    # compute potential of initial distributions, index slice to return 1D vector
    # (note: all columsn are identical in this object by construction)
    phi = eval(sim_params['compute_electric_potential_phi_handle'][x.str])(
        fe, fi, x, vx, 0, sim_params)[:,0]

    # write 1D vector phi = phi(x) to data file
    DECSKS.lib.write.tofiles(
        sim_params['sim_name'],
        n = 0, # time step
        phi = phi
        ) # all cols are identical

# --------------------------------------------------------------------------- #
# *************************************************************************** #
# (4) SIMULATION EXECUTION
# --------------------------------------------------------------------------- #
print 'simulation has started, status updates broadcasted every timestep'
for n in t.stepnumbers:
    fe, fi = DECSKS.lib.split.scheme(
        fe, fi,
        t, x, vx, ax,
        n,
        sim_params
        )

    # write fe, fi to files (./etc/outputs/sim_name/), if phi is computed
    # it is written to file inside lib.fieldsolvers orchestrator
    DECSKS.lib.write.tofiles(sim_name, n, fe = fe, fi = fi)
    DECSKS.lib.status.check_and_clean(t, n, tic)
    toc = time.clock() # stop global clock
# --------------------------------------------------------------------------- #
# *************************************************************************** #
# (5) PRINT COMPLETION NOTICE AND SIMULATION TIME (CLEANUP REPORT TO FOLLOW)
# --------------------------------------------------------------------------- #
simtime = toc - tic
print "simulation completed in %g seconds = %g minutes = %g hours " % (
    simtime,
    simtime/60.,
    simtime/3600.)
# --------------------------------------------------------------------------- #
#=============================================================================#
# END
