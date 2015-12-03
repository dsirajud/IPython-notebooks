import numpy as np
import pylab

def general_maxwellian(f, v, v0 = 0., var = 2., tau = 1., mu = 1.):
    for i in range(f.shape[0]):
        for j in range(f.shape[1]):
            f[i,j] = 1 / np.sqrt(var*np.pi * tau / mu) * np.exp(-(v[j] - v0) ** 2 / (var * tau / mu))
    return f

def maxwellian_exp(f, v, v0 = 0., var = 2., tau = 1., mu = 1.):
    for i in range(f.shape[0]):
        for j in range(f.shape[1]):
            f[i,j] = np.exp(-(v[j] - v0) ** 2 / (var * tau / mu))
    return f

def bump_on_tail(f,x,v):
    for i in range(f.shape[0]):
        for j in range(f.shape[1]):
            f[i,j] = 1 / np.sqrt(2*np.pi) * (1 + 0.04 * np.cos(0.3*x[i]))

    f1 = np.zeros_like(f)
    f2 = np.zeros_like(f)

    f *= (0.9 * maxwellian_exp(f1, v, v0 = 0, var = 2.) +  \
          0.2 * maxwellian_exp(f2, v, v0 = 4.5, var = 1/4.))

    return f

def ion_density(f, x, A):
    f_x = np.zeros_like(f)
    for i in range(f.shape[0]):
        for j in range(f.shape[1]):
            f_x[i,j] = A*x[i] ** 2

    f *= f_x

    return f

ax, bx = -2*np.pi / .3, 2*np.pi / .3
Nx = 512
x = np.linspace(ax, bx, Nx)
xwidth = (bx - ax) / (Nx - 1)

av, bv = -8.0, 8.0
Nv = 512
v = np.linspace(av, bv, Nv)
vwidth = (bv - av) / (Nv - 1)

X, V = np.meshgrid(x, v)

tau = 1/30.
mu = 1836.15267389
v0 = -4.5

fe = np.zeros([Nx, Nv])
fe = bump_on_tail(fe, x, v)

fi = np.zeros([Nx, Nv])
fi_max = np.zeros_like(fi)
fi_max = general_maxwellian(fi_max, v, v0 = 0, var = 10000., tau = tau, mu = mu)
#fi = ion_density(fi_max,x,A)

#pylab.pcolormesh(X, V, fe.T)
pylab.pcolormesh(X,V, fi_max.T)
pylab.clim(0,0.38) # for Landau test case
pylab.colorbar()
pylab.axis([-10, 10, -8, 8])
pylab.show()
pylab.figure()
