import numpy as np

total_dims = [3,4]
active_dims = [3,4]

ax = -2 * np.pi / .3
bx = 2 * np.pi / .3
avx = -8.0
bvx =  8.0
N = 5



def assemble_spectral_derivative_operator(total_dims, active_dims,
                                          ax, bx, avx, bvx,
                                          N):
    """Returns a dictionary Xi with key/value pairs:

        Xi['x'] -- (ndarray, ndim=3, dtype=complex)
        Xi['vx'] -- (ndarray, ndim=3, dtype=complex)

    Each of these matrices correspond to a matrix with entries

      $$Xi = ((Xi)_{q,i,j}) = 1j * (Delta z xi_{i,j})^q$$

    USAGE NOTE: computing Xi * Ff, where Ff is numpy.fft.fft(f)
    and f.shape = (x.N, vx.N) produces the Fourier transform
    of the derivative coefficients $F[d] equiv D$, where
    D[q,i,j] corresponds to the qth order derivative coefficient
    at a phase space location [i,j]. The method
    lib.derivatives.trigonometric3D takes the row-wise inverse
    transform so that the tensor d[q,i,j] is generated.
    """

    Xi = {}
    xi = {}

    Lx = float(bx - ax)
    xwidth = float((bx - ax) / (total_dims[0] - 1))

    Lvx = float(bvx - avx)
    vxwidth = float((bvx - avx) / (total_dims[1] - 1))

    # build x-advection objects
    N_xi, N_cols = active_dims[0], active_dims[1]
    wave_index = np.arange(N_xi)
    xi_x = np.where(wave_index <= N_xi / 2,
              2*np.pi*wave_index / Lx,
              2*np.pi*(wave_index - N_xi) / Lx)

    xi['x'] = xi_x
    xi_2D = np.outer(xi_x, np.ones(N_cols)) # wavenumbers enumerated by row, copies in Nv columns

    dn = np.arange(1,N).reshape(N-1,1,1)
    Xi['x'] = (1j * xwidth * xi_2D) ** dn

    # build v-advection objects
    N_xi, N_rows = active_dims[1], active_dims[0]
    wave_index = np.arange(N_xi)
    xi_vx = np.where(wave_index <= N_xi / 2,
              2*np.pi*wave_index / Lvx,
              2*np.pi*(wave_index - N_xi) / Lvx)

    xi['vx'] = xi_vx
    xi_2D = np.outer(np.ones(N_rows), xi_vx) # wavenumbers enumerated by row, copies in Nv columns

    xi_2D = np.transpose (xi_2D) # dimensions are (vx.N, x.N) now to match v-advection cases
    dn = np.arange(1,N).reshape(N-1,1,1)
    Xi['vx'] = (1j * vxwidth * xi_2D) ** dn

    return Xi, xi
