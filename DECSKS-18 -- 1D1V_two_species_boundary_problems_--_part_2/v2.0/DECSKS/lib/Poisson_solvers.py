import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt

def Poisson_6th_LDBC_LNBC(n, Nx = 24, a = 0., b = 1., LDBC = 0., LNBC = 1.):
    """6th order LTE finite difference Poisson solver for two Neumann conditions
    for a suitably well-posed problem. The conditions are enforced as Cauchy
    data at the left boundary, where the second Neumann condition on the right
    hand side will be satisfied by consequence provided the problem meets the
    Neumann solubility condition on the Poisson equation

    inputs:
    n -- (ndarray, ndim = 1) density vector
    Nx -- (int) number of grid points
    a -- (float) left boundary of x-domain
    b -- (float) right boundary of x-domain

    outputs:
    phi -- (ndarray, ndim = 1) numerical solution to potential
    """

    L = float(b - a)
    dx = L / (Nx - 1)

    # augmented density vector N = (0, n), shape = (N+1,)

    N = np.zeros(Nx)

    # assemble differencing matrix

    D = np.zeros((Nx,Nx))

    # LDBC row, (row 0)
    D[0,0] = 1.

    # LNBC row, (row 1)
    D[1,0] = -97/10.
    D[1,1] = 16.
    D[1,2] = -10
    D[1,3] = 5.
    D[1,4] = -3/2.
    D[1,5] = 1/5.

    # Poisson's equation rows
    for i in range(2,Nx):
        D[i,i-2] = 1
        D[i,i-1] = -2
        D[i,i] = 1


    # Assemble FD coefficient matrix on n: B
    B = np.zeros((Nx,Nx))
    for i in range(1,B.shape[0]):
        if i == 1:
            B[i,i-1] = 317 / 240.
            B[i,i] = -133/120.
            B[i,i+1] = 187 / 120.
            B[i,i+2] = -23 / 20.
            B[i,i+3] = 109 / 240.
            B[i,i+4] = -3/40.

        if i == 2:
            B[i, i-2] = 3 / 40.
            B[i, i-1] = 209 / 240.
            B[i,i] = 1 / 60.
            B[i,i+1] = 7 / 120.
            B[i,i+2] = -1 / 40.
            B[i,i+3] = 1 / 240.

        elif 3 <= i <= Nx-2:
            B[i,i-3] = -1/240.
            B[i,i-2] = 1/10.
            B[i,i-1] = 97/120.
            B[i,i] = 1/10.
            B[i,i+1] = -1/240.

        elif i == Nx-1:
            B[i,i-5] = 1/240.
            B[i,i-4] = -1/40.
            B[i,i-3] = 7/120.
            B[i,i-2] = 1/60.
            B[i,i-1] = 209/240.
            B[i,i] = 3/40.


    phi_BC = np.zeros(Nx)

    phi_BC[0] = LDBC
    phi_BC[1] = 6 * dx * LNBC
    b = dx ** 2 * B.dot(n) + phi_BC
    phi = LA.solve(D,b)

    return phi

import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt

def Poisson_6th_RDBC_RNBC(n, Nx = 24, a = 0., b = 1., RDBC = 0., RNBC = 1.):
    """6th order LTE finite difference Poisson solver for two Neumann conditions
    for a suitably well-posed problem. The conditions are enforced as Cauchy
    data at the left boundary, where the second Neumann condition on the right
    hand side will be satisfied by consequence provided the problem meets the
    Neumann solubility condition on the Poisson equation

    inputs:
    n -- (ndarray, ndim = 1) density vector
    Nx -- (int) number of grid points
    a -- (float) left boundary of x-domain
    b -- (float) right boundary of x-domain

    outputs:
    phi -- (ndarray, ndim = 1) numerical solution to potential
    """

    L = float(b - a)
    dx = L / (Nx - 1)

    # augmented density vector N = (0, n), shape = (N+1,)

    N = np.zeros(Nx)

    # assemble differencing matrix

    D = np.zeros((Nx,Nx))

    # LDBC row, (row 0)
    D[-1,-1] = 1.

    # LNBC row, (row 1)
    D[-2,-1] = -97/10.
    D[-2,-2] = 16.
    D[-2,-3] = -10
    D[-2,-4] = 5.
    D[-2,-5] = -3/2.
    D[-2,-6] = 1/5.

    # Poisson's equation rows
    for i in range(Nx-2):
        D[i,i] = 1
        D[i,i+1] = -2
        D[i,i+2] = 1


    # Assemble FD coefficient matrix on n: B
    B = np.zeros((Nx,Nx))
    for i in range(B.shape[0]):
        if i == 0:
            B[i,i] = 3/40.
            B[i,i+1] = 209/240.
            B[i,i+2] = 1/60.
            B[i,i+3] = 7/120.
            B[i,i+4] = -1/40.
            B[i,i+5] = 1/240.

        if 1 <= i < Nx-3:
            B[i,i-1] = -1/240.
            B[i,i] = 1/10.
            B[i,i+1] = 97/120.
            B[i,i+2] = 1/10.
            B[i,i+3] = -1/240.

        elif i == Nx-3:
            B[i,i-3] = 1/240.
            B[i,i-2] = -1/40.
            B[i,i-1] = 7/120.
            B[i,i] = 1/60.
            B[i,i+1] = 209/240.
            B[i,i+2] = 3/40.

        elif i == Nx-2:
            B[i,i+1] = 317 / 240.
            B[i,i] = -133/120.
            B[i,i-1] = 187 / 120.
            B[i,i-2] = -23 / 20.
            B[i,i-3] = 109 / 240.
            B[i,i-4] = -3/40.

        # elif i == Nx - 1: all zeros


    phi_BC = np.zeros(Nx)

    phi_BC[-1] = RDBC
    phi_BC[-2] = -6 * dx * RNBC
    b = dx ** 2 * B.dot(n) + phi_BC
    phi = LA.solve(D,b)

    return phi

import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt

def Poisson_6th_LDBC_RNBC(n, Nx = 24, a = 0., b = 1., LBC = 0., RBC = 1.):
    """6th order LTE finite difference Poisson solver for Dirichlet BCs

    inputs:
    n -- (ndarray, ndim = 1) density vector
    Nx -- (int) number of grid points
    a -- (float) left boundary of x-domain
    b -- (float) right boundary of x-domain

    outputs:
    phi -- (ndarray, ndim = 1) numerical solution to potential
    """

    L = float(b - a)
    dx = L / (Nx - 1)

    # augmented density vector N = (0, n), shape = (N+1,)

    N = np.zeros(Nx)

    # assemble differencing matrix

    D = np.zeros((Nx,Nx))

    # RDBC row
    D[0,0] = 1.

    # LNBC row
    D[-1,-1] = -97/10.
    D[-1,-2] = 16.
    D[-1,-3] = -10
    D[-1,-4] = 5.
    D[-1,-5] = -3/2.
    D[-1,-6] = 1/5.

    # Poisson's equation rows
    for i in range(1,Nx-1):
        D[i,i-1] = 1
        D[i,i] = -2
        D[i,i+1] = 1

    # Assemble FD coefficient matrix on n: B
    B = np.zeros((Nx,Nx))
    for i in range(B.shape[0]):
        # i == 0 row contains all zeros

        if i == 1:

            B[i, i-1] = 3 / 40.
            B[i, i] = 209 / 240.
            B[i,i+1] = 1 / 60.
            B[i,i+2] = 7 / 120.
            B[i,i+3] = -1 / 40.
            B[i,i+4] = 1 / 240.

        elif 2 <= i <= Nx-3:

            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == Nx-2:

            B[i,i+1] = 3 / 40.
            B[i,i] = 209 / 240.
            B[i,i-1] = 1 / 60.
            B[i,i-2] = 7 / 120.
            B[i,i-3] = -1 / 40.
            B[i,i-4] = 1 / 240.

        if i == Nx-1:
            B[i,i-5] = -3/40.
            B[i,i-4] = 109 / 240.
            B[i,i-3] = -23 / 20.
            B[i,i-2] = 187 / 120.
            B[i,i-1] = -133/120.
            B[i,i] = 317 / 240.

    phi_BC = np.zeros(Nx)

    phi_BC[0] = LBC
    phi_BC[-1] = -6 * dx * RBC
    b = dx ** 2 * B.dot(n) + phi_BC
    phi = LA.solve(D,b)

    return phi

def Poisson_6th_LNBC_RDBC(n, Nx = 24, a = 0., b = 1., LBC = 0., RBC = 1.):
    """6th order LTE finite difference Poisson solver for Dirichlet BCs

    inputs:
    n -- (ndarray, ndim = 1) density vector
    Nx -- (int) number of grid points
    a -- (float) left boundary of x-domain
    b -- (float) right boundary of x-domain

    outputs:
    x -- (ndarray, ndim = 1) mesh points (abscissa)
    phi -- (ndarray, ndim = 1) numerical solution to potential
    phi_exact -- (ndarray, ndim = 1) exact solution (chosen)
    """

    L = float(b - a)
    dx = L / (Nx - 1)

    # augmented density vector N = (0, n), shape = (N+1,)

    N = np.zeros(Nx)

    # assemble differencing matrix

    D = np.zeros((Nx,Nx))

    # LNBC row
    D[0,0] = -97/10.
    D[0,1] = 16.
    D[0,2] = -10
    D[0,3] = 5.
    D[0,4] = -3/2.
    D[0,5] = 1/5.

    # RDBC row
    D[-1,-1] = 1.

    # Poisson's equation rows
    for i in range(1,Nx-1):
        D[i,i-1] = 1
        D[i,i] = -2
        D[i,i+1] = 1


    # Assemble FD coefficient matrix on n: B
    B = np.zeros((Nx,Nx))
    for i in range(B.shape[0]):
        if i == 0:
            B[i,i] = 317 / 240.
            B[i,i+1] = -133/120.
            B[i,i+2] = 187 / 120.
            B[i,i+3] = -23 / 20.
            B[i,i+4] = 109 / 240.
            B[i,i+5] = -3/40.

        elif i == 1:

            B[i, i-1] = 3 / 40.
            B[i, i] = 209 / 240.
            B[i,i+1] = 1 / 60.
            B[i,i+2] = 7 / 120.
            B[i,i+3] = -1 / 40.
            B[i,i+4] = 1 / 240.

        elif 2 <= i <= Nx-3:

            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == Nx-2:

            B[i,i+1] = 3 / 40.
            B[i,i] = 209 / 240.
            B[i,i-1] = 1 / 60.
            B[i,i-2] = 7 / 120.
            B[i,i-3] = -1 / 40.
            B[i,i-4] = 1 / 240.

        # else i == Nx-1: row of zeros

    phi_BC = np.zeros(Nx)

    phi_BC[0] = 6 * dx * LBC
    phi_BC[-1] = RBC

    b = dx ** 2 * B.dot(n) + phi_BC
    phi = LA.solve(D,b)

    return phi

def Poisson_6th_LDBC_RDBC(n, Nx = 24, a = 0, b = 1, LBC = 0, RBC = 0):
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
    # domain parameters
    dx = float(b - a) / (Nx - 1)
    
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



def Poisson_6th_PBC(n, Nx = 24, a = 0, b = 1, LBC = None, RBC = None):
    """6th order LTE finite difference Poisson solver for Dirichlet BCs

    inputs:
    n -- (ndarray, ndim=1), len(n) = Nx-1, i.e. size active domain
    Nx -- (int) number of grid points
    a -- (float) left boundary of x-domain
    b -- (float) right boundary of x-domain

    outputs:
    error -- (float) L2 norm
    x -- (ndarray, ndim = 1) mesh points (abscissa)
    phi -- (ndarray, ndim = 1) numerical solution to potential
    phi_exact -- (ndarray, ndim = 1) exact solution (chosen)
    """
    # domain parameters
    L = float(b - a)
    dx = L / (Nx - 1)

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = B*n + phi_BC

    # Assemble FD coefficient matrix on phi: D
    D = np.zeros([Nx-1,Nx-1])
    for i in range(Nx-1):
        if i == 0:
            D[i,i] = -2
            D[i,i+1] = 1
            D[i,Nx-2] = 1

        elif i == Nx-2:
            D[i,i] = -2
            D[i,i-1] = 1
            D[i,0] = 1
        else:
            D[i,i-1] = 1
            D[i,i] = -2
            D[i,i+1] = 1

    # Assemble FD coefficient matrix on n: B
    B = np.zeros([Nx-1,Nx-1])
    for i in range(Nx-1):
        if i == 0:
            B[i,Nx-3] = -1/240.
            B[i,Nx-2] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == 1:
            B[i,Nx-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif 1 < i < Nx-3:
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == Nx-3:
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,0] = -1/240.

        elif i == Nx-2:
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,0] = 1/10.
            B[i,1] = -1/240.

    # label the RHS a b = dx ** 2 * B*n, here n[Nx-1] = n[0] not included so len(n) = Nx-1
    b = dx ** 2 * B.dot(n)

    # solve D*phi = b
    phi = LA.solve(D,b)
    phi_total = np.zeros(Nx)
    phi_total[:Nx-1] = phi

    # periodic BC
    phi_total[Nx-1] = phi[0]

    # PBCs do not produce unique solutions but a family of solutions with arbitrary integration constant
    # C = phi_avg over active region i = 0, 1, ... , Nx - 2

    phi_avg = np.sum(phi_total[:Nx-1])*dx / L

    phi_total -= phi_avg

    return phi_total

