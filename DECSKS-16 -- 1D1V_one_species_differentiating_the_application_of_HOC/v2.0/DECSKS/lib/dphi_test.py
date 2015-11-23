import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA
from dphi import assemble_finite_difference_weight_matrix_single_derivative


a, b = 0, 1
L = float(b - a)

grids = np.array([12, 24, 48, 96, 192, 384, 768, 1536])
num_grids = len(grids)
error_norm = np.zeros(num_grids)


print "Nx          error          order"
for i in range(num_grids):

    x = np.linspace(a,b,grids[i])
    dx = float(b - a) / (grids[i] - 1)
    phi = np.sin(2*np.pi*x)
    dphi_exact = 2*np.pi * np.cos(2*np.pi*x)
    W = assemble_finite_difference_weight_matrix_single_derivative(grids[i], 1, 6)
    W /= dx ** 1
    dphidx = W.dot(phi)

    error_norm[i] = LA.norm(dphi_exact - dphidx, 2) * np.sqrt(dx / L)

    if i == 0:
        print "%d          %2.5e          ----" % (grids[i], error_norm[i])
    else:
        order = np.log2(error_norm[i-1] / error_norm[i])
        print "%d          %2.5e          %g" % (grids[i], error_norm[i], order)

plt.plot(x, dphi_exact, label = 'exact', linewidth = 2)
plt.hold('on')
plt.plot(x, dphidx, '-o', label = 'numerical', linewidth = 2)
plt.grid()
plt.legend(loc = 'best')
plt.show()
