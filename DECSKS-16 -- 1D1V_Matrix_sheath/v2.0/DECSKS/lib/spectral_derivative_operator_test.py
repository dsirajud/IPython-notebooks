import numpy as np

total_dims = [101, 3]

active_dims = [100, 2]
ax = 0
bx = np.pi

avx = 0
bvx = 0
N = 4

def assemble_spectral_derivative_operator(total_dims, active_dims,
                                          ax, bx, avx, bvx,
                                          N):

    L = float(bx - ax)
    width = float((bx - ax) / (total_dims[0] - 1))

    # TODO x-differentiation, haven't figured out what v-differentaition looks like yet
    N_xi, N_cols = active_dims[0], active_dims[1]

    # for v-differentiation, number of wavenumbers N_xi = active_dims[1], N_rows (or
    # N_cols if we do transpose operations) = active_dims[0]
    
    wave_index = np.arange(N_xi)
    xi = np.where(wave_index <= N_xi / 2,
              2*np.pi*wave_index / L,
              2*np.pi*(wave_index - N_xi) / L)

    xi = np.outer(xi, np.ones(N_cols)) # wavenumbers enumerated by row, copies in Nv columns

    dn = np.arange(1,N).reshape(N-1,1,1)
    return (1j * width * xi) ** dn
