import matplotlib.pyplot as plt
import numpy as np
from convergence_routines import *

from convergence_routines import FD_derivative_matrix_formulation, domain, df3

x, dx, L = domain(_Nx = 2688)
L2error, df3_approx = FD_derivative_matrix_formulation(_dn = 3, _p = 3, _Nx = 2688)
df3_exact = df3(x)

plt.plot(x,df3_exact, label = 'exact df3', linewidth = 3)
plt.hold('on')
plt.plot(x,df3_approx, label = 'approx df3', linewidth = 1, color = "red")


# compare with the function whose derivative this is

df2_exact = df2(x)
plt.plot(x,df2_exact * np.abs(np.min(df3_approx)) / np.abs(np.min(df2_exact)), label = 'exact df2', linewidth = 1, color = "cyan")

plt.hold('off')
plt.legend(loc = 'best')
plt.grid()
plt.show()
