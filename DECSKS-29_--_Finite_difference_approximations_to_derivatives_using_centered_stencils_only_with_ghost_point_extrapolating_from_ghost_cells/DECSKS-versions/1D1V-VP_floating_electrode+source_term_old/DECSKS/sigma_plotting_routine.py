
# inside timestep loop

    sim_params['sigma_n']['x']['lower'][n] = sim_params['sigma']['x']['lower']
    sim_params['sigma_n']['x']['upper'][n] = sim_params['sigma']['x']['upper']

# after timestep loop, final plot summary

sigma_n_left = sim_params['sigma_n']['x']['lower']
sigma_n_right = sim_params['sigma_n']['x']['upper']

trange = np.arange(t.N+1)
plt.plot(trange,sigma_n_left, linewidth = 2, label = r'$\sigma (t, x= -10)$')
plt.plot(trange,sigma_n_right, linewidth = 2, label = r'$\sigma (t, x= +10)$')
plt.grid()
plt.xlabel(r'time step $n$', fontsize = 18)
plt.ylabel(r'$\sigma (t,x)$', fontsize = 18)
plt.legend(loc = 'best')
plt.figure()


phi_left = 1/2. * sim_params['sigma_n']['x']['lower'] # E = -1/2 sigma, phi = 1/2 sigma, here sigma = ni - ne
phi_right = 1/2. * sim_params['sigma_n']['x']['upper']
plt.plot(trange,phi_left, linewidth = 2, label = r'$\phi (t, x= -10)$')
plt.plot(trange,phi_right, linewidth = 2, label = r'$\phi (t, x= +10)$')
plt.grid()
plt.xlabel(r'time step $n$', fontsize = 18)
plt.ylabel(r'$\phi (t,x)$', fontsize = 18)
plt.legend(loc = 'best')

plt.show()
