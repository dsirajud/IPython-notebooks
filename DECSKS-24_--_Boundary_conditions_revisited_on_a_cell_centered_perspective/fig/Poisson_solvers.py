import numpy as np
import numpy.linalg as LA

def LDBC_RDBC(n, Nx = 24, dx = .5, LBC = 0, RBC = 0):
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

    # Dirichlet boundary conditions (BC) vector
    phi_DBC = np.zeros(Nx)
    phi_DBC[0] = -2 * LBC
    phi_DBC[-1] = -2 * RBC

    # label the RHS a b = dx ** 2 * B*n + phi_DBC
    b = dx ** 2 * B.dot(n) + phi_DBC

    # solve D*phi = b
    phi = LA.solve(D,b)

    return phi

def LNBC_RDBC(n, Nx = 24, dx = .5, LBC = 0, RBC = 0):
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
            D[i,i] = -1
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

    # boundary conditions (BC) vector
    phi_BC = np.zeros(Nx)
    phi_BC[0] = dx * LBC
    phi_BC[-1] = -2 * RBC

    # label the RHS a b = dx ** 2 * B*n + phi_DBC


    # solve D*phi = b
    phi = LA.solve(D,b)

    return phi

def PBC(n, ax = 0., bx = 2., Nx = 24, dx = .5, LBC = 0, RBC = 0):
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
            D[i,-1] = 1
            D[i,i] = -2
            D[i,i+1] = 1
        elif i == Nx-1:
            D[i,0] = 1
            D[i,i] = -2
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

    # boundary conditions (BC) vector
    phi_BC = np.zeros(Nx)
    phi_BC[0] = 0
    phi_BC[-1] = 0

    # label the RHS a b = dx ** 2 * B*n + phi_DBC
    b = dx ** 2 * B.dot(n) + phi_BC

    # solve D*phi = b
    phi = LA.solve(D,b)
    phi_avg = np.sum(phi)*dx / (bx - ax)
    phi -= phi_avg

    return phi

def PBCv2(n, ax = 0., bx = 2., Nx = 24, dx = .5, LBC = 0, RBC = 0):
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
            D[i,-1] = 1
            D[i,i] = -2
            D[i,i+1] = 1
        elif i == Nx-1:
            D[i,0] = 1
            D[i,i] = -2
            D[i,i-1] = 1
        else:
            D[i,i-1] = 1
            D[i,i] = -2
            D[i,i+1] = 1

    # Assemble FD coefficient matrix on n: B
    B = np.zeros([Nx,Nx])
    for i in range(Nx):
        if i == 0:
            B[i,-2] = -1/240.
            B[i,-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        if i == 1:
            B[i,-1] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.


        elif 1 < i < Nx-2:
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == Nx-2:
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,0] = -1/240.

        elif i == Nx-1:
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,0] = 1/10.
            B[i,1] = -1/240.

    # boundary conditions (BC) vector
    phi_BC = np.zeros(Nx)
    phi_BC[0] = 0
    phi_BC[-1] = 0

    # label the RHS a b = dx ** 2 * B*n + phi_DBC
    b = dx ** 2 * B.dot(n) + phi_BC

    # solve D*phi = b
    phi = LA.solve(D,b)
    phi_avg = np.sum(phi)*dx / (bx - ax)
    phi -= phi_avg

    return phi


