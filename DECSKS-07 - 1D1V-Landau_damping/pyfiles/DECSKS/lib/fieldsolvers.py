import numpy as np
import numpy.linalg as LA
import DECSKS

def Gauss(ni, f, x, v, n):
    """Computes self-consistent electric field E by solving Poisson's equation
    using FFT/IFFT.

    inputs:
    ni -- (float) uniform background density of ions,
                  in the future can take an input fi, to compute ni
    f -- (ndarray, dim=2) electron density fe(x,v,n) at time step t^n
                          used to compute ne(x,n) at time step t^n
    x -- (instance) spatial variable
    v -- (instance) velocity variable
    n -- (int) time step number, t^n


    outputs:
    E -- (ndarray,dim=1) electric field, E(x) at time t^n
    """
    ne = DECSKS.lib.density.single_integration(f[n,:x.N,:v.N], of = x, wrt = v)

    xi    = np.zeros(x.N)
    E_hat = np.zeros(x.N, dtype = complex) # container


    # define wave indices
    for r in range(x.N):
        if r <= x.N/2 :
            xi[r] = 2*np.pi*r / x.L
        else:
            xi[r] = 2*np.pi*(r - x.N) / x.L

    # total charge density, n(x), len(n) = Nx

    n = ni - ne

    N = np.fft.fft(n)
    A    = max(N)
    eps   = 2.0e-15
    xi_min = A*eps
    for r in range(x.N):
        if np.abs(N[r]) < xi_min:
            N[r] = 0

    # E_hat[0] = 0 from periodic BCs, i.e. because N[0] = 0, quasineutrality
    # equivalently, E_hat[0] is the DC field, of which there is none in
    # a quasineutral system of charged particles, only flucutations are present
    # E_hat[0] = 0 already from E_hat vector initialization

    for each in range(1,len(xi)):

        E_hat[each] = 1 / (1j*xi[each]) * N[each]         # Electric field in Fourier space


    E = np.real(np.fft.ifft(E_hat))    # Electric field in configurational space


    return E

def Poisson(ni, f, x, v, n):
    """6th order LTE finite difference Poisson solver for Dirichlet BCs

    inputs:
    ni -- (float) uniform background density of ions,
                  in the future can take an input fi, to compute ni
    f -- (ndarray, dim=2) electron density fe(x,v,n) at time step t^n
                          used to compute ne(x,n) at time step t^n
    x -- (instance) spatial variable
    v -- (instance) velocity variable
    n -- (int) time step number, t^n

    outputs:
    phi -- (ndarray,dim=1) electric potential, phi(x) at time t^n
    E -- (ndarray,dim=1) electric field, E(x) at time t^n

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
    # pass full array f at t = t^n, boundaries included
    ne = DECSKS.lib.density.single_integration(f[n,:,:], of = x, wrt = v)
    n = ne - ni

    # the boundary values from phi_exact are seen to be:
    phi_a = 0
    phi_b = np.sin(np.pi**2)

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = B*n + phi_BC

    # Assemble FD coefficient matrix on phi: D
    D = np.zeros([Nx,Nx])
    for i in range(x.Ngridpoints):
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




# ........................................................................... #
