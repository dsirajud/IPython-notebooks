import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA
import Poisson_solvers

if __name__ == '__main__':

    Nx = 12
    a,b = 0, 2
    L = b - a
    dx = float(L) / Nx

    x = np.zeros(Nx)
    x[0] = a + dx / 2.
    for i in range(Nx):
        x[i] = x[0] + i*dx

    x_exact = np.linspace(a, b, 200)
    x_exact_points = np.linspace(a,b,Nx)

    phi_exact = - 1 / np.pi ** 2 * np.sin(np.pi * x_exact)

    x_edges = np.zeros(Nx+1)

    x_edges[0] = 0
    x_edges[1:-1] = x[1:] - dx/2
    x_edges[-1] = 2.

    phi_exact_edges = -1 / np.pi ** 2 * np.sin(np.pi * x_edges)

    n = np.sin(np.pi * x)
    phi = Poisson_solvers.PBC(n, ax = 0., bx = 2., Nx = Nx, dx = dx, LBC = None, RBC = None)

    plt.figure()
    plt.plot(x_exact, phi_exact, '--', markersize = 25, lw = 2,label = 'exact', color = 'cornflowerblue')
    plt.plot(x_edges, phi_exact_edges, '|', markersize = 25, lw = 3,label = 'exact values at cell edges', color = 'b')
    plt.plot(x, phi, 'o', lw = 2, label = 'PBC solution')
    #    plt.annotate(r'$\partial_x\phi_{exact} (0) = 1$', xy=(0, phi[0]), xytext=(0.2, phi[0]+0.5), ha = 'center',
    #            arrowprops=dict(facecolor='black', shrink=0.05),
    #            )

    #    plt.annotate(r'$\phi (%g) = 0$' % x[0], xy=(x[0], phi[0]), xytext=(x[0], 0.8), ha = 'center',
    #            arrowprops=dict(facecolor='black', shrink=0.05),
    #            )


    #    plt.annotate(r'$\phi_{exact} (2) = 3$' % x[-1], xy=(x_exact[-1],phi_exact[-1]), xytext=(x_exact[-1]-.1, 1.6), ha = 'right',
    #            arrowprops=dict(facecolor='black', shrink=0.05),
    #            )


    plt.grid()
    #    plt.axis([-.2, 2.1, -.2, 3.2])
    plt.legend(loc = 'best',)
    plt.savefig('Poisson_PBC_solver_perf.png')
