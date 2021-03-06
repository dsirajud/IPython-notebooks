class phasespace_var:
    def __init__(self, Nz, az, bz, name = 'x'):
        self.N = Nz
        self.az = az
        self.bz = bz
        self.L = float(bz - az)
        self.width = self.L / (self.N - 1)
        self.gridvalues = np.linspace(az, bz, Nz)

        self.str = name

import numpy as np

# constants
epsi = 0.04
k = 0.3
vT1 = 2.
mu1 = 0.
vT2 = 1/4.
mu2 = 4.5
A = 1 / np.sqrt(2*np.pi)
B = 0.9
C = 0.2

N = 21

Nx, Nv = 768, 1536
ax, bx = -2*np.pi / k, 2*np.pi / k
av, bv = -8., 8.

x = phasespace_var(Nx, ax, bx, name = 'x')
v = phasespace_var(Nv, av, bv, name = 'v')
X,V = np.meshgrid(x.gridvalues,v.gridvalues)

f = np.zeros([x.N,v.N])
f = A*(1 + epsi*np.cos(k*X))*(B*np.exp( -(V - mu1) ** 2 / vT1) + C*np.exp( -(V - mu2) ** 2 / vT2))

f = f.T

 
