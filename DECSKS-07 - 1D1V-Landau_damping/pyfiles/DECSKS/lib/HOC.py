import numpy as np
from scipy.misc import factorial

def beta(p,a):
    """Returns the pth coefficient beta.

    inputs:
    p -- (int) pth coefficient
    a -- (float) uncorrected fractional distance of MC

    output:

    beta_p -- pth beta coefficient for given vel. or accel.
    """

    beta_p = a ** (p+1) / factorial(p+1)

    if p > 0:

        for q in range(p):

            if a >= 0:

                beta_p -= beta(q,a) / factorial(p + 1 - q)

            else: # a < 0

                beta_p -= (-1) ** (p+q) * beta(q,a) / factorial(p + 1 - q)

    return beta_p
# ........................................................................... #
def beta_m(a, B, N):
    A = np.zeros([N,N])
    alpha = np.zeros([N,1])
    for i in range(N):
        alpha[i,0] = a ** (i+1) / factorial(i+1)
        for j in range(i+1):
            A[i,j] = B[i-j]/factorial(i-j)
    return A.dot(alpha)
# ........................................................................... #
def kernel(z):
    """Returns a Gaussian enveloped window filtered wave number
    in Fourier space.

    inputs:
    z -- (instance) phase space variable

    output:
    K -- filtered wave number K(xi)
    """
    sigma = 4
    mu    = 0.0

    K = np.sinc(z.cells /z.width) * np.exp(- (z.cells - mu)**2 / (2*sigma**2))

    return K
# ........................................................................... #
