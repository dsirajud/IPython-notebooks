import matplotlib.pyplot as plt
import numpy as np
from convergence_routines import *

Nx = 2688
x, dx, L = domain(_Nx = Nx)
L2error, df5_approx = FD_derivative_matrix_formulation(_dn = 5, _p = 3, _Nx = Nx)
df5_exact = df5(x)

plt.plot(x,df5_exact, label = 'exact df5', linewidth = 3)
plt.hold('on')
plt.plot(x,df5_approx, label = 'approx df5', linewidth = 1, color = "red")


# compare with the function whose derivative this is

df4_exact = df4(x)
plt.plot(x,df4_exact * np.abs(np.min(df5_approx)) / np.abs(np.min(df4_exact)), label = 'exact df4', linewidth = 1, color = "cyan")

plt.hold('off')
plt.legend(loc = 'best')
plt.xlim([0.131, 0.146])
plt.ylim([-3.2e6, -2.4e6])
plt.grid()
plt.show()
