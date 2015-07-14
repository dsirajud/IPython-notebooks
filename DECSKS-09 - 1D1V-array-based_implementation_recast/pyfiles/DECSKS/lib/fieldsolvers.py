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

def Poisson_PBC_6th(ni, f,
                x, v, n):
    """6th order LTE finite difference Poisson solver for periodic BCs

    inputs:
    ni -- (float) uniform background density of ions,
                  in the future can take an input fi, to compute ni
    f -- (ndarray, dim=2) electron density fe(x,v,n) at time step t^n
                          used to compute ne(x,n) at time step t^n
    x -- (instance) spatial variable
    v -- (instance) velocity variable
    n -- (int) time step number, t^n


    outputs:
    phi -- (ndarray,dim=1) scalar potential, phi(x) at time t^n,
           for i = 0, 1, ... , x.N - 1, one full period
    """

    # charge densities
    ne = DECSKS.lib.density.single_integration(f[n,:x.N,:v.N], of = x, wrt = v)
    n = ne - ni

    # form the tensor objects involved in the numerical solution
    #
    #     d^2 phi = n --> D*phi = B*n + phi_BC

    # Assemble FD coefficient matrix on phi: D
    D = np.zeros([x.N,x.N])
    for i in range(x.N):
        if i == 0:         # first row
            D[i,i] = -2
            D[i,i+1] = 1
            D[i,-1] = 1

        elif i == x.N - 1: # last row
            D[i,i] = -2
            D[i,i-1] = 1
            D[i,0] = 1
        else:              # interior rows
            D[i,i-1] = 1
            D[i,i] = -2
            D[i,i+1] = 1

    # Assemble FD matrix B
    B = np.zeros([x.N,x.N])
    for i in range(x.N):
        if i == 0:             # first row
            B[i,-2] = -1/240.
            B[i,-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == 1:           # second row
            B[i,-1] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif 1 < i < (x.N - 2): # 2 <= row < = third before last
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,i+2] = -1/240.

        elif i == (x.N - 2): # second before last row
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,i+1] = 1/10.
            B[i,0] = -1/240.

        elif i == (x.N - 1): # last row
            B[i,i-2] = -1/240.
            B[i,i-1] = 1/10.
            B[i,i] = 97/120.
            B[i,0] = 1/10.
            B[i,1] = -1/240.

    # label the RHS a b = dx ** 2 * B*n
    b = x.width ** 2 * B.dot(n)

    # solve D*phi = b
    phi = LA.solve(D,b)

    # PBCs do not produce unique solutions but a family of solutions with arbitrary integration constant
    # that corresponds to a DC offset phi_avg, recenter so that phi_avg = 0

    phi_avg = np.sum(phi) * x.width / x.L
    phi -= phi_avg

    return phi
# ........................................................................... #
