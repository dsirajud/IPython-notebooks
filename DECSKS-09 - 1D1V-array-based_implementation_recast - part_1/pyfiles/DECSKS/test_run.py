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
sim_params['m_0'] = np.sum(DECSKS.lib.domain.extract_active_grid(f[0,:,:], x, sim_params))

# Current case: uniform background (cold) density of ions,
sim_params['ni'] = DECSKS.lib.density.cold_background(f,x,vx,sim_params)


f_initial = f[0,:,:]

#############################################
# X ADVECTION TEST

# (0) INITIALIZE FINAL DENSITY CONTAINER AND EXTRACT EVOLVED GRID
f_final = np.zeros_like(f_initial)
f_initial = DECSKS.lib.domain.extract_active_grid(f_initial, x, sim_params)

z = x
vz = vx

z.CFL.compute_numbers(z, vz, 2.)

# advection step
z.postpointmesh = DECSKS.lib.convect.advection_step(z)

z.postpointmesh = DECSKS.lib.boundaryconditions.periodic(z.postpointmesh, z.N)



# remap_step
f_old = f_initial # needed for some function call references
c = DECSKS.lib.HOC.correctors(sim_params, z, vz)
d = np.zeros([sim_params['N'], z.N, vz.N])
# zeroeth derivative coefficient is the density itself
d[0,:,:] = f_old
if sim_params['N'] > 1:
    d[1:,:,:] = DECSKS.lib.derivatives.trigonometric3D(f_old, z, sim_params)

# evaluate derivatives q = 0, 1, 2, ... N-1 (zeroeth is density itself)

# calls lib.derivatives.fd or lib.derivatives.fourier based on the
# HOC specified in etc/params.dat, sim_params['derivative_method']
# contains the related function handle as a string
d = eval(sim_params['derivative_method'])(f_old, z, vz, sim_params)

# compute high order fluxes
Uf = DECSKS.lib.convect.flux(
    sim_params,
    f_old,
    z, vz
    )

    # enforce flux limiter to ensure positivity and restrict numerical overflow
Uf = DECSKS.lib.convect.flux_limiter(f_old, Uf, z)

f_k1 = DECSKS.lib.convect.remap_assignment(
    f_old,
    Uf,
    z.postpointmesh[0,:,:],
    z,
    vz,
    index = 'nearest'    # remaps to nearest neighbor index
    )

# remap to contiguous cell center
f_k2 = DECSKS.lib.convect.remap_assignment(
    f_old,
    Uf,
    z.postpointmesh[1,:,:],
    z,
    vz,
    index = 'contiguous'    # remaps to contiguous neighbor of above
    )


f_remapped = f_k1 + f_k2

n = 1 # current time step, just manually putting in number, is actually passed to function

# global check on density conservation
DECSKS.lib.density.global_conservation_check(sim_params, f_remapped, n)

f_final = DECSKS.lib.convect.finalize_density(sim_params, f_remapped, f_final, z)

#############################################
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
sim_params['m_0'] = np.sum(DECSKS.lib.domain.extract_active_grid(f[0,:,:], x, sim_params))

# Current case: uniform background (cold) density of ions,
sim_params['ni'] = DECSKS.lib.density.cold_background(f,x,vx,sim_params)


f_initial = f[0,:,:]
#############################################
# V ADVECTION test

Ex = DECSKS.lib.fieldsolvers.Gauss1D1V(sim_params['ni'], f, x, vx, 0, sim_params)
ax = DECSKS.lib.domain.Setup(sim_params, 'a', 'x')
ax.prepointvaluemesh = -Ex

z = vx
vz = ax

f_final = np.zeros(f_initial.shape)
f_initial = DECSKS.lib.domain.extract_active_grid(f_initial, z, sim_params)

# (1) PUSH MOVING CELLS an integer number of cells

z.CFL.compute_numbers(z, vz, 2.)

# check if velocity being advected, if so then transpose several objects
if z.str[0] == 'v':
    f_initial, z, vz = DECSKS.lib.domain.velocity_advection_prep(f_final, f_initial, z, vz)

# advection step
z.postpointmesh = DECSKS.lib.convect.advection_step(z)
z.postpointmesh = DECSKS.lib.boundaryconditions.periodic(z.postpointmesh, z.N)



# remap_step
f_old = f_initial # needed for some function call references
c = DECSKS.lib.HOC.correctors(sim_params, z, vz)
d = np.zeros([sim_params['N'], z.N, vz.N])
# zeroeth derivative coefficient is the density itself
d[0,:,:] = f_old
if sim_params['N'] > 1:
    d[1:,:,:] = DECSKS.lib.derivatives.trigonometric3D(f_old, z, sim_params)

# evaluate derivatives q = 0, 1, 2, ... N-1 (zeroeth is density itself)

# calls lib.derivatives.fd or lib.derivatives.fourier based on the
# HOC specified in etc/params.dat, sim_params['derivative_method']
# contains the related function handle as a string
d = eval(sim_params['derivative_method'])(f_old, z, vz, sim_params)

# compute high order fluxes
Uf = DECSKS.lib.convect.flux(
    sim_params,
    f_old,
    z, vz
    )

    # enforce flux limiter to ensure positivity and restrict numerical overflow
Uf = DECSKS.lib.convect.flux_limiter(f_old, Uf, z)

f_k1 = DECSKS.lib.convect.remap_assignment(
    f_old,
    Uf,
    z.postpointmesh[0,:,:],
    z,
    vz,
    index = 'nearest'    # remaps to nearest neighbor index
    )

# remap to contiguous cell center
f_k2 = DECSKS.lib.convect.remap_assignment(
    f_old,
    Uf,
    z.postpointmesh[1,:,:],
    z,
    vz,
    index = 'contiguous'    # remaps to contiguous neighbor of above
    )


f_remapped = f_k1 + f_k2

n = 1 # current time step, just manually putting in number, is actually passed to function

# global check on density conservation
DECSKS.lib.density.global_conservation_check(sim_params, f_remapped, n)

f_final = DECSKS.lib.convect.finalize_density(sim_params, f_remapped, f_final, z)



