import _mypath     # adds relative path to sys.path for flexible deployment
import DECSKS
import time
# =========================================================================== #

#rm_plots = int(raw_input('remove ALL plot files after simulation is done (1 = yes, 0 = no)?: '))
#tic = time.clock()

sim_params = DECSKS.lib.read.inputfile('./etc/params.dat')
sim_params['BC'] = 'periodic'
x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
vx = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')
f = DECSKS.lib.density.setup(sim_params, t, x, vx)    # f = f(x_i, v_j, t^n)

# Current case: uniform background (cold) density of ions,
sim_params['ni'] = DECSKS.lib.density.cold_background(f,x,vx,sim_params)

if sim_params['HOC'] == 'FD':
    sim_params['W'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrices(sim_params,x,vx)
    sim_params['W_dn1'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrix_single_derivative(sim_params,x,dn = 1, LTE = 6)

DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,0,x,vx)

#Plot = DECSKS.lib.plots.PlotSetup(f, 0, t, x, v, sim_params)
#Plot(n = 0)

print 'simulation has started, status updates are broadcasted after each timestep'

for n in t.stepnumbers:

    f = DECSKS.lib.split.scheme(
        f,
        t,x,vx,
        n,
        sim_params
        )

    #    Plot = DECSKS.lib.plots.PlotSetup(f, n, t, x, v, sim_params)
    #    Plot(n)

    # calcs performed and outputs written only if "record outputs? = yes"
    # in ./etc/params.dat
    DECSKS.lib.diagnostics.calcs_and_writeout(sim_params,f,n,x,vx)
    DECSKS.lib.status.check_and_clean(t, n, tic, rm_plots)
