import numpy as np
import DECSKS

def Poisson(ni, f, x, v, n):
    """Computes self-consistent electric field E by solving Poisson's equation
    using FFT/IFFT.

    inputs:
    n0 -- (float) uniform background density of ions,
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

# ........................................................................... #
