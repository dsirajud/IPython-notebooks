import matplotlib.pyplot as plt
import numpy as np
from convergence_routines import *

Nx = 2488
x, dx, L = domain(_Nx = Nx)
L2error, df9_approx = FD_derivative_matrix_formulation(_dn = 9, _p = 3, _Nx = Nx)
df9_exact = df9(x)

plt.plot(x,df9_exact, label = 'exact df9', linewidth = 3)
plt.hold('on')
plt.plot(x,df9_approx, label = 'approx df9', linewidth = 1, color = "red")


# compare with the function whose derivative this is

df8_exact = df8(x)
plt.plot(x,df8_exact * np.abs(np.min(df9_approx)) / np.abs(np.min(df8_exact)), label = 'exact df4', linewidth = 1, color = "cyan")

plt.hold('off')
plt.legend(loc = 'best')
plt.grid()
plt.show()
