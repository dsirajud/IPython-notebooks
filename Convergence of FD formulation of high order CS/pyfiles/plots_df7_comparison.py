import matplotlib.pyplot as plt
import numpy as np
from convergence_routines import *

Nx = 1344
x, dx, L = domain(_Nx = Nx)
L2error, df7_approx = FD_derivative_matrix_formulation(_dn = 7, _p = 3, _Nx = Nx)
df7_exact = df7(x)

plt.plot(x,df7_exact, label = 'exact df7', linewidth = 3)
plt.hold('on')
plt.plot(x,df7_approx, label = 'approx df7', linewidth = 1, color = "red")


# compare with the function whose derivative this is

df6_exact = df6(x)
plt.plot(x,df6_exact * np.abs(np.min(df7_approx)) / np.abs(np.min(df6_exact)), label = 'exact df4', linewidth = 1, color = "cyan")

plt.hold('off')
plt.legend(loc = 'best')
plt.grid()
plt.show()
