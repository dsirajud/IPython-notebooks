import numpy as np
import numpy.linalg as LA

def main(Nx = 24, a = 0, b = 1, return_all_results = 'no'):
    """6th order LTE finite difference Poisson solver for Dirichlet BCs

    inputs:
    Nx -- (int) number of grid points
    a -- (float) left boundary of x-domain
    b -- (float) right boundary of x-domain
    return_all_results -- (str) if 'yes', function returns ( x, phi, phi_exact, error)
                          elif 'no', function returns (error)
                          
    outputs:
    error -- (float) L2 norm
    x -- (ndarray, ndim = 1) mesh points (abscissa)
    phi -- (ndarray, ndim = 1) numerical solution to potential
    phi_exact -- (ndarray, ndim = 1) exact solution (chosen)
    """
    # mesh setup
    a = float(a)
    b = float(b)
    L = b - a
    x = np.linspace(a, b, Nx) # spacing dx = (b - a ) / (Nx - 1)
    dx = L / (Nx - 1)

    
    # phi_exact is chosen, Poisson's equation gives n
    n = -np.pi**4 * np.sin(np.pi**2 * x)
    phi_exact = np.sin(np.pi**2 * x)
    
    # the boundary values from phi_exact are seen to be:
    phi_a = 0
    phi_b = np.sin(np.pi**2)
    
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
        if i == 0:
            B[i,i] = 0

        elif i == 1:
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

    # Dirichlet boundary conditions (DBC) vector
    phi_DBC = np.zeros(Nx)
    phi_DBC[0] = phi_a
    phi_DBC[-1] = phi_b

    # label the RHS a b = dx ** 2 * B*n + phi_DBC
    b = dx ** 2 * B.dot(n) + phi_DBC

    # solve D*phi = b
    phi = LA.solve(D,b)
    
    # define normalized root-mean square error, i.e. normalized L2 norm of the error
    error_norm = LA.norm(phi_exact - phi,2) * np.sqrt(dx / L)

    if return_all_results.lower() != 'no':
        return error_norm, x, phi, phi_exact
    else:
        return error_norm
