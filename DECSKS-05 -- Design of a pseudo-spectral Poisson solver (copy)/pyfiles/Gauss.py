import numpy as np

def Gauss(rho, x):
    """Computes self-consistent electric field E by solving Gauss' law
    using FFT/IFFT.

    inputs:
    rho -- (ndarray, ndim = 1)  charge density vector,
    x -- (ndarrya, ndim = 1) configurational variable

    outputs:
    E -- (ndarray,dim=1) electric field, E(x) at time t^n
    """
    Nx = len(x)
    L = np.max(x) - np.min(x)
    dx = L / (Nx - 1)
    
    # create wave number vector in same order as return of numpy.fft.fft
    xi = np.zeros(Nx)
    xi[0:Nx/2 + 1] = 2 * np.pi * np.arange(Nx/2 + 1) / L
    xi[Nx/2 + 1:] = 2 * np.pi * (np.arange(Nx/2 + 1, Nx) - Nx) / L

    # Fluctuating electric field in Fourier space, FE_tilde[0] = 0 already
    FE_tilde = np.zeros(Nx, dtype = complex)
    Frho_tilde = np.fft.fft(rho)
    FE_tilde[1:Nx] = Frho_tilde[1:Nx] / (1j * xi[1:Nx])
    
    # Electric field in configurational space
    E_tilde = np.real(np.fft.ifft(FE_tilde))
    return E_tilde
