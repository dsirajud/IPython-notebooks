import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA

def Poisson_6th_LDBC_RDBC(n, Nx = 24, dx = .5, a = 0, b = 1, LBC = 0, RBC = 0):
    """6th order LTE finite difference Poisson solver for Dirichlet BCs
       over the domain x in [a,b]

    inputs:
    n -- (ndarray, ndim=1) density in Poisson's equation
    Nx -- (int) number of grid points
    phi_a -- (float) potential at left endpoint, x = a
    phi_b -- (float) potential at right endpoint, x = b

    outputs:
    phi -- (ndarray, ndim = 1) numerical solution to potential
    """

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = B*n + phi_BC

    # Assemble FD coefficient matrix on phi: D
    D = np.zeros([Nx,Nx])
    for i in range(Nx):
        if i == 0 or i == Nx-1:
            D[i,i] = 1
        else:
            D[i,i-1] = 1
            D[i,i] = -2
            D[i,i+1] = 1

    # Assemble FD coefficient matrix on n: B
    B = np.zeros([Nx,Nx])
    for i in range(Nx):
        # if i == 0: all zeros (to make row DBC equation independent of density n)

        if i == 1:
            B[i,i-1] = 3/40.
            B[i,i] = 209/240.
            B[i,i+1] = 1/60.
            B[i,i+2] = 7/120.
            B[i,i+3] = -1/40.
            B[i,i+4] = 1/240.

        elif i == Nx-1:
            B[i,i] = 0

        elif 1 < i < Nx-2:
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == Nx-2:
            B[i,i-4] = 1/240.
            B[i,i-3] = -1/40.
            B[i,i-2] = 7/120.
            B[i,i-1] = 1/60.
            B[i,i] = 209/240.
            B[i,i+1] = 3/40.

        #elif i == Nx-1: all zeros (to make row DBC equation independent of density n)

    # Dirichlet boundary conditions (DBC) vector
    phi_DBC = np.zeros(Nx)
    phi_DBC[0] = LBC
    phi_DBC[-1] = RBC

    # label the RHS a b = dx ** 2 * B*n + phi_DBC
    b = dx ** 2 * B.dot(n) + phi_DBC

    # solve D*phi = b
    phi = LA.solve(D,b)

    return phi

def Poisson_6th_LDBC_RDBC_new(n, Nx = 24, dx = .5, a = 0, b = 1, LBC = 0, RBC = 0):
    """6th order LTE finite difference Poisson solver for Dirichlet BCs
       over the domain x in [a,b]

    inputs:
    n -- (ndarray, ndim=1) density in Poisson's equation
    Nx -- (int) number of grid points
    phi_a -- (float) potential at left endpoint, x = a
    phi_b -- (float) potential at right endpoint, x = b

    outputs:
    phi -- (ndarray, ndim = 1) numerical solution to potential
    """

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = B*n + phi_BC

    # Assemble FD coefficient matrix on phi: D
    D = np.zeros([Nx,Nx])
    for i in range(Nx):
        if i == 0:
            D[i,i] = -3
            D[i,i+1] = 1
        elif i == Nx-1:
            D[i,i] = -3
            D[i,i-1] = 1
        else:
            D[i,i-1] = 1
            D[i,i] = -2
            D[i,i+1] = 1

    # Assemble FD coefficient matrix on n: B
    B = np.zeros([Nx,Nx])
    for i in range(Nx):
        if i == 0: 
            B[i,i] = 317/240.
            B[i,i+1] = -133/120.
            B[i,i+2] = 187/120.
            B[i,i+3] = -23/20.
            B[i,i+4] = 109/240.
            B[i,i+5] = -3/40.

        if i == 1:
            B[i,i-1] = 3/40.
            B[i,i] = 209/240.
            B[i,i+1] = 1/60.
            B[i,i+2] = 7/120.
            B[i,i+3] = -1/40.
            B[i,i+4] = 1/240.

        elif 1 < i < Nx-2:
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == Nx-2:
            B[i,i-4] = 1/240.
            B[i,i-3] = -1/40.
            B[i,i-2] = 7/120.
            B[i,i-1] = 1/60.
            B[i,i] = 209/240.
            B[i,i+1] = 3/40.

        elif i == Nx-1:
            B[i,i] = 317/240.
            B[i,i-1] = -133/120.
            B[i,i-2] = 187/120.
            B[i,i-3] = -23/20.
            B[i,i-4] = 109/240.
            B[i,i-5] = -3/40.


    # Dirichlet boundary conditions (DBC) vector
    phi_DBC = np.zeros(Nx)
    phi_DBC[0] = -2 * LBC
    phi_DBC[-1] = -2 * RBC

    # label the RHS a b = dx ** 2 * B*n + phi_DBC
    b = dx ** 2 * B.dot(n) + phi_DBC

    # solve D*phi = b
    phi = LA.solve(D,b)

    return phi

if __name__ == '__main__':

    print "first we reuse the LDBC/RDBC solver with the proper interpretation as the gridpoints as cell centers"

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

    phi_exact = 1/2. * x_exact ** 2 + 1/2. * x_exact

    x_edges = np.zeros(Nx+1)

    x_edges[0] = 0
    x_edges[1:-1] = x[1:] - dx/2
    x_edges[-1] = 2.

    phi_exact_edges = 1/2. * x_edges ** 2 + 1/2. * x_edges

    n = np.ones(Nx)
    phi = Poisson_6th_LDBC_RDBC(n, Nx, dx, a, b, LBC = 0., RBC = 3.)


    plt.figure()
    plt.plot(x_exact, phi_exact, '--', markersize = 25, lw = 2,label = 'exact', color = 'cornflowerblue')
    plt.plot(x_edges, phi_exact_edges, '|', markersize = 25, lw = 3,label = 'exact values at cell edges', color = 'b')
    plt.plot(x, phi, 'o', lw = 2, label = 'old LDBC/RDBC solution')
    plt.annotate(r'$\phi_{exact} (0) = 0$', xy=(0, 0), xytext=(0, 0.5), ha = 'center',
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

    plt.annotate(r'$\phi (%g) = 0$' % x[0], xy=(x[0], phi[0]), xytext=(x[0], 0.8), ha = 'center',
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

    plt.annotate(r'$\phi (%g) = 3$' % x[-1], xy=(x[-1],phi[-1]), xytext=(1.2, phi[-1]),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

    plt.annotate(r'$\phi_{exact} (2) = 3$' % x[-1], xy=(x_exact[-1],phi_exact[-1]), xytext=(x_exact[-1], 1.5), ha = 'right',
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

    plt.grid()
    plt.axis([-.2, 2.1, -.2, 3.2])
    plt.legend(loc = 'best',)
    plt.savefig('Old_Poisson_LDBC_RDBC_solver_perf.png')

    print "plot saved as %s" % 'Old_Poisson_LDBC_RDBC_solver_perf.png'
    print ""
    print "next, we try the new solution"


    n = np.ones(Nx)

    phi_new = Poisson_6th_LDBC_RDBC_new(n, Nx, dx, a, b, LBC = 0., RBC = 3.)

    plt.figure()
    plt.plot(x_exact, phi_exact, '--', markersize = 25, lw = 2,label = 'exact', color = 'cornflowerblue')
    plt.plot(x_edges, phi_exact_edges, '|', markersize = 25, lw = 3,label = 'exact values at cell edges', color = 'b')
    plt.plot(x, phi, 'o', lw = 2, label = 'old LDBC/RDBC solution')
    plt.plot(x, phi_new, 'D', lw = 2, label = 'new LDBC/RDBC solution')
    plt.annotate(r'$\phi (0) = 0$', xy=(0, 0), xytext=(0, 0.5), ha = 'center',
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
    plt.annotate(r'$\phi (%g) = 3$' % x[-1], xy=(x[-1],phi[-1]), xytext=(1.2, phi[-1]),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

    plt.annotate(r'$\phi (2) = 3$' % x[-1], xy=(x_exact[-1],phi_exact[-1]), xytext=(x_exact[-1], 1.5), ha = 'right',
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

    plt.grid()
    plt.axis([-.2, 2.1, -.2, 3.2])
    plt.legend(loc = 'best',)
    plt.savefig('New_Poisson_LDBC_RDBC_solver_perf.png')

    print "plot saved as %s" % 'New_Poisson_LDBC_RDBC_solver_perf.png'
    print ""
