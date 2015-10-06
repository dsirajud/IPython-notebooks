# run this in ipython from DECSKS/
import _mypath     # adds relative path to sys.path for flexible deployment
import DECSKS
import numpy as np
# =========================================================================== #

#rm_plots = int(raw_input('remove ALL plot files after simulation is done (1 = yes, 0 = no)?: '))
#tic = time.clock()

sim_params = DECSKS.lib.read.inputfile('./etc/params.dat')

x = DECSKS.lib.domain.Setup(sim_params, var = 'x')
vx = DECSKS.lib.domain.Setup(sim_params, var = 'v', dim = 'x')
t = DECSKS.lib.domain.Setup(sim_params, var = 't')
f = DECSKS.lib.density.setup(sim_params, t, x, vx)    # f = f(x_i, v_j, t^n)

# store total mass for conservation checks
sim_params['m_0'] = np.sum(DECSKS.lib.convect.extract_active_grid(sim_params, f[0,:,:]))

# Current case: uniform background (cold) density of ions,
sim_params['ni'] = DECSKS.lib.density.cold_background(f,x,vx,sim_params)


f_initial = f[0,:,:]
Ex = DECSKS.lib.fieldsolvers.Gauss(sim_params['ni'], f, x, vx, 0)
ax = DECSKS.lib.domain.Setup(sim_params, 'a', 'x')
ax.prepointvaluemesh = -Ex

# (0) INITIALIZE FINAL DENSITY CONTAINER AND EXTRACT EVOLVED GRID
f_final = np.zeros_like(f_initial)
f_initial = DECSKS.lib.convect.extract_active_grid(sim_params, f_initial)

z = x
vz = vx

z.CFL.compute_numbers(z, vz, 2.)
z, k_x = DECSKS.lib.convect.advection_step(f_initial, z, vz)
z.postpointmesh = DECSKS.lib.boundaryconditions.periodic(z.postpointmesh, z.N)
k_x = DECSKS.lib.boundaryconditions.periodic(k_x, z.N)



# (1) PUSH MOVING CELLS an integer number of cells
z = vx
vz = ax

z.CFL.compute_numbers(z, vz, 200.)
z, k_vx = DECSKS.lib.convect.advection_step(f_initial, z, vz)
z.postpointmesh = DECSKS.lib.boundaryconditions.periodic(z.postpointmesh, z.N)
k_vx = DECSKS.lib.boundaryconditions.periodic(k_vx, z.N)
k_vx_T = np.transpose(k_vx, axes = (0,2,1)) # flips each 2D matrix at a time

f_old = f_initial


# analyze x advection for now

z = x
vz = vx

# The following W and W_dn1 excerpt is from main.py

if sim_params['HOC'] == 'FD':
    sim_params['W'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrices(sim_params,x,vx)
    sim_params['W_dn1'] = DECSKS.lib.derivatives.assemble_finite_difference_weight_matrix_single_derivative(sim_params,x,dn = 1, LTE = 6)

Uf = DECSKS.lib.convect.flux(
    sim_params,
    f_old,
    z, vz
    )
