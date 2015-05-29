import matplotlib.pyplot as plt
import numpy as np
from convergence_routines import *

Nx = 2688
x, dx, L = domain(_Nx = Nx)
L2error, df8_approx = FD_derivative_matrix_formulation(_dn = 8, _p = 3, _Nx = Nx)
df8_exact = df8(x)

plt.plot(x,df8_exact, label = 'exact df8', linewidth = 3)
plt.hold('on')
plt.plot(x,df8_approx, label = 'approx df8', linewidth = 1, color = "red")


# compare with the function whose derivative this is

df7_exact = df7(x)
plt.plot(x,df7_exact * np.abs(np.min(df8_approx)) / np.abs(np.min(df7_exact)), label = 'exact df4', linewidth = 1, color = "cyan")

plt.hold('off')
plt.legend(loc = 'best')
plt.grid()
plt.show()
