# after Plot(n = 0) call

#create and store densities at a chosen x, here we choose x = 0, or x.gridpoints[240 / 2] = x.gridpoints[120]
fe_v = np.zeros([t.N+1, fe.shape[2]])
fi_v = np.zeros([t.N+1, fi.shape[2]])

fe_v[0,:] = fe[0,120,:]
fi_v[0,:] = fi[0,120,:]

#inside simulation loop:
    fe_v[n,:] = fe[n,120,:]
    fi_v[n,:] = fi[n,120,:]


# after simulation completes, plot all

for n in range(t.N+1):
    plt.plot(vx.gridvalues,fe_v[n,:], linewidth = 2, label = r'$f_e(t = %g,x=0,v_x)$' % t.times[n])
    plt.grid()
    plt.xlabel(r'$v_x$', fontsize = 18)
    plt.ylabel(r'$f_i(x)$', fontsize = 18)
    plt.legend(loc = 'best')

plt.figure()
for n in range(t.N+1):
    plt.plot(vx.gridvalues,fi_v[n,:], linewidth = 2, label = r'$f_i(t = %g,0,v_x)$' % t.times[n])
    plt.grid()
    plt.xlabel(r'$v_x$', fontsize = 18)
    plt.ylabel(r'$f_i(x)$', fontsize = 18)
    plt.legend(loc = 'best')

plt.show()
