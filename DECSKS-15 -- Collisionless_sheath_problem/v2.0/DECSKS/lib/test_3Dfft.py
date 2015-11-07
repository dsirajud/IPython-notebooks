import numpy as np
import matplotlib.pyplot as plt

def trigonometric2D(Y,
                    Nx,
                    xwidth,
                    xi
                    ):
    """Computes derivatives of density for derivative coeffs d in FN methods

    inputs:
    f -- (ndarray, dim=1) f(z1,z2=const,z3=const,..., t = n-1)
    z -- (instance) phase space variable being convected
    q -- (int) order of desired derivative
    xi -- (ndarray, dim=2) 1D wave numbers stacked in columns for column-wise
          fourier transform.
    sim_params -- (dict) simulation parameters
    K -- (ndarray, ndim=1) Windowed Fourier transform kernel

    outputs:
    d -- (ndarray, dim=1) qth derivative coefficient for all z.prepoints
    """

    Ff = np.fft.fft(Y)

    #    eps = np.finfo(np.float)()

    #    Ff = np.where(Ff < eps, 0, Ff)

    # calc first derivative only
    D = (1j*xi*xwidth) * Ff
    d = np.real(np.fft.ifft(D))

    return d

# tier 1
x = np.linspace(0,2*np.pi, 100)
y = np.sin(x)
y2 = np.sin(2*x)
y3 = np.sin(3*x)

# tier 2
yy = np.cos(x)
yy2 = np.cos(2*x)
yy3 = np.cos(3*x)


Y = np.zeros([2, x.shape[0], 3])

for n in range(1,3+1):
    Y[0,:,n-1] = np.sin(n*x)
    Y[1,:,n-1] = np.cos(n*x)


# calculate wave numbers

Nx = len(x)
Nv = Y.shape[2]

xwidth = x[1] - x[0]
L = x[-1] - x[0]

wave_index = np.arange(Nx)

xi = np.where(wave_index <= Nx / 2,
              2*np.pi*wave_index / L,
              2*np.pi*(wave_index - Nx) / L)

xi = np.outer(xi, np.ones(Nv))

d = trigonometric2D(Y,
                    Nx,
                    xwidth,
                    xi)

