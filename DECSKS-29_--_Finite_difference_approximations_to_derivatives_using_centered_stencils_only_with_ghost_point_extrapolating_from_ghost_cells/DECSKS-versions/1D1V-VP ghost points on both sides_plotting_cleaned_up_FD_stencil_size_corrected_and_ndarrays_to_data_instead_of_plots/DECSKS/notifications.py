#======================================================================#
# DECSKS.notifications                                                 #
#----------------------------------------------------------------------#
# prints to screen various parameters to be used in the simulation     #
# this was designed just to make sure the input file is being read in  #
# consistently                                                         #
#----------------------------------------------------------------------#

import numpy as np

def announce(fe, fi, x, vx, sim_params, print_to_screen = True):

    if print_to_screen:
        print "===================================================================================================="
        print "The following information has been requested by the user in DECSKS.notification.announce to be printed "
        print "----------------------------------------------------------------------------------------------------\n"

        print "DECSKS will use the following Poisson solver: %s" % sim_params['compute_electric_potential_phi_handle'][x.str]
        print "\n"

        print "boundary conditions specified include:\n"

        print "BC type on f on x:  %s" %  sim_params['BC']['f']['x']['type']
        print "BC type on f on vx: %s" %  sim_params['BC']['f']['vx']['type']
        print "\n"
        print "BC (left)  on f(x):  %s" %  sim_params['BC']['f']['x']['lower']
        print "BC (right) on f(x):  %s" %  sim_params['BC']['f']['x']['upper']
        print '\n'
        print "BC (left)  on f(vx): %s" %  sim_params['BC']['f']['vx']['lower']
        print "BC (right) on f(vx): %s" %  sim_params['BC']['f']['vx']['upper']

        IQ = np.sum(-fe + fi) * vx.width * x.width
        print "\ninitial charge in system IQ = %g\n" % IQ
        print "====================================================================================================\n\n"

    return None
