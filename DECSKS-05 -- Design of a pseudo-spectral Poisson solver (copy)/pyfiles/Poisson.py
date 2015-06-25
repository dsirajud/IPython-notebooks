import numpy as np

def Poisson(rho, x, return_phi = 'yes', return_dphi = 'no'):
    """Computes self-consistent potential phi and/or dphi by 
    solving Poisson's eq.using FFT/IFFT.

    inputs:
    rho -- (ndarray, ndim = 1)  charge density vector,
    x -- (ndarrya, ndim = 1) configurational variable

    outputs:
    phi and/or dphi -- (ndarray,dim=1) electric field, E(x) at time t^n
    """
    Nx = len(x)
    L = np.max(x) - np.min(x)
    dx = L / (Nx - 1)
    
    # create wave number vector in same order as return of numpy.fft.fft
    xi = np.zeros(Nx)
    xi[0:Nx/2 + 1] = 2 * np.pi * np.arange(Nx/2 + 1) / L
    xi[Nx/2 + 1:] = 2 * np.pi * (np.arange(Nx/2 + 1, Nx) - Nx) / L

    # Fluctuating dphi  in Fourier space, Fdphi_tilde[0] = 0 already
    Frho_tilde = np.fft.fft(rho)
    
    if return_dphi.lower() == 'yes':
        Fdphi_tilde = np.zeros(Nx, dtype = complex)
        Fdphi_tilde[1:Nx] = -Frho_tilde[1:Nx] / (1j * xi[1:Nx])
        dphi_tilde = np.real(np.fft.ifft(Fdphi_tilde))
    
    if return_phi.lower() == 'yes':
        Fphi_tilde = np.zeros(Nx, dtype = complex)
        Fphi_tilde[1:Nx] =  Frho_tilde[1:Nx] / (xi[1:Nx]) ** 2
        phi_tilde = np.real(np.fft.ifft(Fphi_tilde))
    
    if return_dphi.lower() == 'yes' and return_phi.lower() == 'yes':
        return dphi_tilde, phi_tilde
    
    elif return_dphi.lower() != 'yes' and return_phi.lower() == 'yes':
        return phi_tilde
    
    elif return_dphi.lower() == 'yes' and return_phi.lower() != 'yes':
        return dphi_tilde
