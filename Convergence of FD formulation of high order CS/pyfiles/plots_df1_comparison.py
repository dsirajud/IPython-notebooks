import matplotlib.pyplot as plt
import numpy as np
from convergence_routines import *

x, dx, L = domain(_Nx = 2688)
L2, df1_approx = FD_derivative_matrix_formulation(_dn = 1, _p = 3, _Nx = 2688)
df1_exact = df1(x)

plt.plot(x,df1_exact, label = 'exact df1', marker = 'o', linewidth = 2)
plt.hold('on')
plt.plot(x,df1_approx, label = 'approx df1', linewidth = 2)
plt.hold('off')

plt.legend(loc = 'best')
plt.grid()
plt.show()
