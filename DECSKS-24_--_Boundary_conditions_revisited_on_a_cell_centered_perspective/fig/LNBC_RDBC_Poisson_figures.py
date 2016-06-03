import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA
import Poisson_solvers

if __name__ == '__main__':

    Nx = 6
    a,b = 0, 2
    L = b - a
    dx = float(L) / Nx

    x = np.zeros(Nx)
    x[0] = a + dx / 2.
    for i in range(Nx):
        x[i] = x[0] + i*dx

    x_exact = np.linspace(a, b, 200)
    x_exact_points = np.linspace(a,b,Nx)

    phi_exact = 1/2. * x_exact ** 2 + x_exact - 1

    x_edges = np.zeros(Nx+1)

    x_edges[0] = 0
    x_edges[1:-1] = x[1:] - dx/2
    x_edges[-1] = 2.

    phi_exact_edges = 1/2. * x_edges ** 2 + x_edges - 1

    n = np.ones(Nx)
    phi = Poisson_solvers.LNBC_RDBC(n, Nx = Nx, dx = dx, LBC = 1, RBC = 3)

    plt.figure()
    plt.plot(x_exact, phi_exact, '--', markersize = 25, lw = 2,label = 'exact', color = 'cornflowerblue')
    plt.plot(x_edges, phi_exact_edges, '|', markersize = 25, lw = 3,label = 'exact values at cell edges', color = 'b')
    plt.plot(x, phi, 'o', lw = 2, label = 'LNBC/RDBC solution')
    plt.annotate(r'$\partial_x\phi_{exact} (0) = 1$', xy=(0, phi[0]), xytext=(0.2, phi[0]+0.5), ha = 'center',
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

    #    plt.annotate(r'$\phi (%g) = 0$' % x[0], xy=(x[0], phi[0]), xytext=(x[0], 0.8), ha = 'center',
    #            arrowprops=dict(facecolor='black', shrink=0.05),
    #            )


    plt.annotate(r'$\phi_{exact} (2) = 3$' % x[-1], xy=(x_exact[-1],phi_exact[-1]), xytext=(x_exact[-1]-.1, 1.6), ha = 'right',
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

    # plot the line with slope dphi(0) = 1, whose y intercept is -1
    LNBC = 1*x_exact - 1
    plt.plot(x_exact, LNBC, '--', lw = 2,label = r'slope $= 1$', color = 'k')
    plt.grid()
    #    plt.axis([-.2, 2.1, -.2, 3.2])
    plt.legend(loc = 'best',)
    plt.savefig('Poisson_LNBC_RDBC_solver_perf.png')
