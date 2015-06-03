import matplotlib.pyplot as plt
import numpy as np
from convergence_routines import *

from convergence_routines import FD_derivative_matrix_formulation, domain, df2

x, dx, L = domain(_Nx = 336) # Nx = 336 shows convergent FD derivative from above output
L2error, df2_approx = FD_derivative_matrix_formulation(_dn = 2, _p = 3, _Nx = 336)
df2_exact = df2(x)

plt.plot(x,df2_exact, label = 'exact df2', marker = 'o', linewidth = 2)
plt.hold('on')
plt.plot(x,df2_approx, label = 'approx df2', linewidth = 2)
plt.hold('off')

plt.legend(loc = 'best')
plt.grid()
plt.show()
